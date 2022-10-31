# %load ./FunctionsDB.py


import copy
from MainFunctions import writepickle, initray
from rdkit import Chem  # Importing RDKit
from FindFunctionalGroups import identify_functional_groups as IFG
import os
import dask.dataframe as dd
import dask.bag as db
import pandas as pd
from typing import List, Dict, Union, Tuple
from AnalgCompds import getCarrierFrags0
import modin.pandas as mpd

#%% Building basic substance database


def info(molfile: str) -> Union[Tuple[Chem.rdchem.Mol, str], str]:
    """
    Processes a molfile address of a species and returns either a tuple of the SMILES and the mol object, or the address if the mol object is not valid.

    Args:
        molfile (str): molfile address of the species.

    Returns:
        Union[Tuple[Chem.rdchem.Mol, str], str]: Either a tuple of the SMILES and the mol object, or the address if the mol object is not valid.
    """
    mol = Chem.MolFromMolFile(molfile)
    if mol:
        Chem.SanitizeMol(mol)
        mol.UpdatePropertyCache(strict=False)
        smiles = Chem.MolToSmiles(mol)
        return mol, smiles
    else:
        return molfile


def basic(ID: str, folder: str) -> Dict:
    """
    Based on the species ID and folder location, returns a dictionary of species information by processing the molfile.
    The molfile needs to be located in folder+os.sep+ID address, and will be processed by RDKit.
    Args:
        ID (str): Reaxys ID of the species.
        folder (str): Folder address of the species mol file on the server

    Returns:
        Dict: Dictionary of species information with SubstanceID, MolFileAddress, and Error if present
    """

    if str.isdecimal(ID):
        molfileaddress = folder + os.sep + ID
        try:
            res = info(molfileaddress)
        except Exception as e:
            error = e
            compaddrs = {
                "SubstanceID": int(ID),
                "MolFileAddress": molfileaddress,
                "Error": error,
            }
        else:
            if type(res) == tuple:
                smiles = res[1]
                compaddrs = {
                    "SubstanceID": int(ID),
                    "MolFileAddress": molfileaddress,
                    "Smiles": smiles,
                }
            else:
                error = "Valence error"
                compaddrs = {
                    "SubstanceID": int(ID),
                    "MolFileAddress": molfileaddress,
                    "Error": error,
                }
    else:
        compaddrs = {}
    return compaddrs


def basicgroup(molfilelist: List[str], folder: str) -> List[Dict]:
    """
    Invokes basic for each ID in molfilelist and returns a list of dictionaries of species information.

    Args:
        molfilelist (List[str]): List of Reaxys IDs to be processed
        folder (str): Folder address

    Returns:
        List[Dict]: List of dictionaries of species information
    """
    return [basic(ID, folder) for ID in molfilelist]


def substancedblist(folderName: str, partitions: int) -> List[Dict]:
    """
    Invoke basicgroup for folder and returns a list of dictionaries of species information.

    Args:
        folderName (str): Folder address
        partitions (int): Number of partitions to split folder contents into.

    Returns:
        List[Dict]: List of dictionaries of species information
    """
    dem = os.listdir(folderName)
    b = db.from_sequence(dem, npartitions=partitions)
    dflist = b.map_partitions(basicgroup, folderName).compute()
    return dflist


#%% Adding fragment smarts column accounting for mixtures
def getfrags(
    series: pd.Series, expand: int = 1, resFormat: str = "smiles"
) -> Union[List[str], str]:  # natoms changed to expand
    """
    Returns carrier fragments of a species given a dataframe row or pandas series and a number of atoms to expand from identified functional groups.
    Accordingly processes mixtures (eg. salts)

    Args:
        series (pd.Series): pandas series of a row of the dataframe
        expand (int, optional): Number of neighboring atoms to expand functional group to. Defaults to 1.
        resFormat (str, optional): Format of the returned results. Defaults to "smiles".

    Returns:
        Union[List[str], str]: If no errors in smiles, returns a tuple of the list of fragment smiles/smarts present in species. Otherwise,
        returns "Error"
    """
    smiles = series["Smiles"]
    if smiles == "Error":
        return "Error"
    if (
        series[">1 Compound"] == True
    ):  # This compound is a mixture. Need to split and apply getcarrierfrags to each smiles
        fraginfo = getmixturefrags(smiles, expand=expand, resFormat=resFormat)
    else:
        try:
            fraginfo = getCarrierFrags0(smiles, expand=expand, resFormat=resFormat)
        except Exception:
            return "Error"
    if type(fraginfo) != list:
        fraginfo = [fraginfo]
    return fraginfo


def createfragdb(
    substancedb: Union[str, pd.DataFrame],
    expand: int = 1,
    ncpus: int = 16,
    restart: bool = True,
    resFormat: str = "smiles",
    explode: bool = True,
    fragdb: pd.DataFrame = None,
) -> pd.DataFrame:
    """
    Generates fragment dataframe from a substance dataframe by applying the getfrags function to each row.

    Args:
        substancedb (Union[str, pd.DataFrame]): Either a path to a pickle file of a pandas dataframe of the substance database, or a pandas dataframe of the substance database.
        expand (int, optional): Number of neighboring atoms to expand functional group to. Defaults to 1.
        ncpus (int, optional): Number of CPUs to use. Defaults to 16.
        restart (bool, optional): Controls if ray cluster already exists, will restart. Defaults to True.
        resFormat (str, optional): Format of the returned results. Defaults to "smiles".
        explode (bool, optional): Controls if the returned dataframe is exploded (on FragmentSmiles or FragmentSmarts). Defaults to True.
        fragdb (pd.DataFrame, optional): If provided, will not run the fragment function but instead will process/cleanup the dataframe. Defaults to None.

    Returns:
        pd.DataFrame: Fragment dataframe with the new column FragmentSmiles or FragmentSmarts
    """
    if resFormat == "smiles":
        columns = "FragmentSmiles"
    elif resFormat == "smarts":
        columns = "FragmentSmarts"
    if fragdb is None:
        if isinstance(substancedb, str):
            substancedb = pd.read_pickle(substancedb)
        if ncpus > 1:
            initray(restart=restart, num_cpus=ncpus)
            substancedbdis = mpd.DataFrame(substancedb)
        else:
            substancedbdis = substancedb
        fragdat = substancedbdis.apply(
            getfrags, axis=1, expand=expand, resFormat=resFormat, result_type="reduce"
        )
        fragdat = pd.Series(data=fragdat.values, index=fragdat.index)
        fragdb = pd.DataFrame(data=fragdat, index=fragdat.index, columns=[columns])
    if fragdb.index.name != "SubstanceID" or fragdb.index.names:
        if fragdb.index.name or fragdb.index.names:
            fragdb.reset_index(inplace=True)
        fragdb.set_index("SubstanceID", inplace=True)
    if substancedb.index.name != "SubstanceID" or substancedb.index.names:
        if substancedb.index.name or substancedb.index.names:
            substancedb.reset_index(inplace=True)
        substancedb.set_index("SubstanceID", inplace=True)
    fragdb["Smiles"] = substancedb.loc[substancedb.index.isin(fragdb.index)]["Smiles"]
    fragdb[">1 Compound"] = substancedb[substancedb.index.isin(fragdb.index)][
        ">1 Compound"
    ]
    # fragdb=fragdb.loc[fragdb[columns]!='Error']
    if explode:
        fragdb = fragdb.explode(columns)
    if "count" not in fragdb.columns:
        fragdb["count"] = fragdb.groupby([fragdb.index, columns])[columns].transform(
            "count"
        )
    fragdb.reset_index(inplace=True)
    fragdb.drop_duplicates(inplace=True)
    fragdb.set_index([columns, "SubstanceID"], inplace=True)
    return fragdb


def createfragfreq(fragdb: Union[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Returns a fragment frequency dataframe

    Args:
        fragdb (Union[str,pd.DataFrame]): Either a path to a pickle file of a pandas dataframe of the fragment database, or a pandas dataframe of the fragment database.

    Returns:
        pd.DataFrame: Fragment frequency dataframe
    """
    if isinstance(fragdb, str):
        fragdb = pd.read_pickle(fragdb)
    if fragdb.index.names or fragdb.index.name:
        fragdb.reset_index(inplace=True)
    freq = (
        fragdb[["SubstanceID", "FragmentSmiles"]]
        .groupby(["FragmentSmiles"])["SubstanceID"]
        .transform("count")
    )
    freqtable = fragdb[["FragmentSmiles"]]
    freqtable["Frequency"] = freq
    freqtable.drop_duplicates(inplace=True)
    freqtable["Popularity"] = (
        freqtable["Frequency"]
        .div(freqtable["Frequency"].sum(axis=0), axis=0)
        .multiply(100)
        .round()
    )
    return freqtable


def getfragpartition(partition, natoms):
    return partition.apply(getfrags, natoms=natoms, axis=1)


#%% Adding fragment smiles column
def getfragsmiles(fragsmarts: Union[str, List[str]]) -> Union[str, List[str]]:
    """
    Converts fragment smarts to smiles

    Args:
        fragsmarts (Union[str,List[str]]): Fragment smarts either string or list of strings

    Returns:
        Union[str,List[str]]: Returns fragment smiles either as a string , list of strings or "Error" if input is "Error"
    """
    if type(fragsmarts) == list:
        fragsmiles = [
            Chem.MolToSmiles(Chem.MolFromSmarts(fragsmart)) for fragsmart in fragsmarts
        ]
    elif fragsmarts == "Error":
        return "Error"
    else:
        fragsmiles = Chem.MolToSmiles(Chem.MolFromSmarts(fragsmarts))
    return fragsmiles


def getsmiles(series: pd.Series) -> Union[str, List[str]]:
    """
    Row version of getfragsmiles

    Args:
        series (pd.Series): Row containing fragment smarts

    Returns:
        Union[str, List[str]]: Returns fragment smiles either as a string , list of strings or "Error" if input is "Error"
    """
    fragsmarts = series["FragmentSmarts"]
    return getfragsmiles(fragsmarts)


def getsmilespartition(partition):
    return partition.apply(getsmiles, axis=1)


#%% Creating mixture column (True if mixture, False if not mixture, Error if smiles not present)


def mixtures(smiles: str) -> Union[bool, str]:
    """
    Returns if a species is a mixture.

    Args:
        smiles (str): SMILES string of the species.

    Returns:
        Union[bool,str]: True if the species is a mixture, False if not, and Error if the species has no SMILES.
    """
    if smiles == "Error":
        return "Error"
    elif len(smiles.split(".")) > 1:
        return True
    else:
        return False


def findMixtures(series: pd.Series) -> Union[bool, str]:
    """
    Takes in a pandas series/row corresponding to a species, and returns if the species is a mixture.

    Args:
        series (pd.Series): Pandas series or row corresponding to a species.

    Returns:
        Union[bool,str]: True if the species is a mixture, False if not, and Error if the species has no SMILES.
    """
    smiles = series["Smiles"]
    return mixtures(smiles)


def findMixturespartition(partition):
    return partition.apply(findMixtures, axis=1)


#%% Changing fragment entries in mixture rows


def getmixturefrags(
    mixsmiles: str, expand: int = 1, resFormat: str = "smarts", addHs: bool = True
) -> Union[str, List[str]]:
    """
    Returns the carrier fragments of a mixture (SMILES has '.' in it) by splitting it into constituent species

    Args:
        mixsmiles (str): Mixture SMILES string
        expand (int, optional): Number of neighboring atoms to expand functional group to. Defaults to 1.
        resFormat (str, optional): Result type, either 'smarts' or 'smiles. Defaults to 'smarts'.
        addHs (bool, optional): If true, hydrogens will be considered in generation of fragment, recommended as true
        otherwise terminal atoms of a molecule are not differentiated with other atoms. Defaults to True.

    Returns:
        Union[str,List[str]]: List of strings of SMARTS or SMILES representing carrier fragments of the mixture. If workflow fails, returns 'Error'.
    """
    try:
        res = []
        for smiles in mixsmiles.split("."):
            reslist = getCarrierFrags0(
                smiles, expand=expand, resFormat=resFormat, addHs=addHs
            )
            if type(reslist) != list:
                res += [reslist]
            else:
                res += reslist
    except Exception:
        return "Error"
    else:
        return res


def getMixturefrags(series: pd.Series, expand: bool = 1) -> Union[str, List[str]]:
    """
    Row version of getmixturefrags

    Args:
        series (pd.Series): Row containing mixture SMILES
        expand (bool, optional): Number of neighboring atoms to expand functional group to. Defaults to 1.

    Returns:
        Union[str,List[str]]: List of strings of SMARTS or SMILES representing carrier fragments of the mixture. If workflow fails, returns 'Error'.
    """
    mixsmiles = series[
        "Smiles"
    ]  # add .values[0] if column is a multiindex, otherwise droplevel = 1 to remove list
    return getmixturefrags(mixsmiles, expand=expand)


def getMixturefragspartition(partition, expand=1):
    return partition.apply(getMixturefrags, expand=expand, axis=1)


#%% Joining columns to a dataframe


def joindf(
    seriesdf: pd.DataFrame, DB: pd.DataFrame, explodeDB: Union[str, Tuple] = None
) -> pd.DataFrame:
    """
    Joins seriesdf (dataframe to add) to DB (dataframe to be added to). Optionally, can explode DB into separate rows for each entry in DB based on
    columns specified in explodeDB.

    Args:
        seriesdf (pd.DataFrame): Dataframe to add.
        DB (pd.DataFrame): Dataframe to be added to.
        explodeDB (Union[str,Tuple], optional): Either a string column or tuple of multiple columns. Defaults to None.

    Returns:
        pd.DataFrame: Dataframe with seriesdf added to DB.
    """
    if seriesdf.index.name != DB.index.name or seriesdf.index.names != DB.index.names:
        if DB.index.name or DB.index.names:
            DB.reset_index(inplace=True)
        if seriesdf.index.name or seriesdf.index.names:
            seriesdf.reset_index(inplace=True)
    DB = DB.join(seriesdf)
    if explodeDB:
        DB = DB.explode(explodeDB)
    return DB


def createspecfreqdb(
    substancedb: Union[str, pd.DataFrame],
) -> pd.DataFrame:
    """
    Creates a dataframe with species and frequency of occurrence in a database.

    Args:
        substancedb (Union[str,pd.DataFrame]): Database to be analyzed. Can be either a filename or a dataframe.

    Returns:
        pd.DataFrame: Dataframe with species and frequency of occurrence in the database.
    """
    if isinstance(substancedb, str):
        substancedb = pd.read_pickle(substancedb)
    species = substancedb.groupby("Smiles").size().reset_index(name="Frequency")
    species.sort_values(by="Frequency", ascending=False, inplace=True)
    species.reset_index(drop=True, inplace=True)
    return species


#%% Deprecated..see DB_Reaxys.ipynb
def buildfragdb(
    sdb=None,
    sdbd=None,
    writesdbd=False,
    sdbdc=None,
    fragseries=None,
    natoms=None,
    writefragseries=False,
    fdbm=None,
    fdb=None,
    dfdb=None,
    fragsmiles=False,
    mixtures=False,
    mixturefrags=False,
    index=None,
    writefdb=False,
    writedfdb=False,
):

    # Note: cluster and client must be initiated for this function to work. Substance database, either dask or pandas
    # should also be loaded. It is advised to avoid loading dask dataframes due to high memory usage. It is recommended to
    # persist and create dask dataframes outside this function, otherwise overheads may be added. Never reset dask dataframe
    # index or else it will reindex all partitions...always start with pandas, reindex and then convert to dask when large computations
    # need to be done.

    # sdb = substance database, sdbd = substance database dask, sdbdc = substance database dask cleaned
    # (remove smile errors), # fragseries = series containing fragment information, natoms =  size of fragment,
    # writefragseries = True if write to file false otherwise, fdbm = fragment database master, pandas version of sdbdc,
    # fdb = final fragment database exploded and unindexed pfdb = final pandas fragment dataframe, multiindexed,
    # dfdb = final dask fragment dataframe

    # Step 1: Create dask dataframe from pandas substance database (output: sdbd)

    if sdb and not type(sdbd) == dd.core.DataFrame:
        sdbd = dd.from_pandas(sdb, npartitions=16)
        sdbd = client.persist(sdbd)
        print("Dask substance dataframe created and persisted")
        if writesdbd:
            sdbd.to_parquet(writesdbd)
            print("Dask substance dataframe written to file: " + writesdbd)
        if not natoms and not fragseries:
            print("Please specify size of fragments that should be retrieved")
            return sdbd

    # Step 2: Clean dask data frame and select only substance ID and smiles column, removing errors (output: sdbdc).

    if sdbd and not type(sdbdc) == dd.core.DataFrame:
        sdbd = client.persist(sdbd)
        sdbdc = sdbd.reset_index()[["SubstanceID", "Smiles", ">1 Compound"]]
        sdbdc = sdbdc[sdbdc.Smiles != "Error"]
        sdbdc = client.persist(sdbdc)
        print("Cleaned dask substance dataframe created and persisted")
        if not natoms and not type(fragseries) == pd.core.series.Series:
            print("Please specify size of fragments that should be retrieved")
            return sdbdc

    # Step 3: Scrape 16 million compounds, and extract series of active fragments for each (output: fragseries)

    if (
        not type(fragseries) == pd.core.series.Series
        and not type(fdb) == pd.core.frame.DataFrame
    ):
        if not type(sdbdc) == dd.core.DataFrame:
            return "Please include a cleaned dask substance dataframe for fragment retrieval"
        if not natoms:
            return "Please specify size of fragments retrieved"
        sdbdc = client.persist(sdbdc)
        #         if natoms==0:
        #             name='ActiveFragmentSmarts'
        #         else:
        #            filename='CarrierFragments'+'(n='+str(natoms)+')'
        name = "FragmentSmarts"
        fragseries = sdbdc.map_partitions(
            getfragpartition, natoms=natoms, meta=(name, "O")
        ).compute()
        print("fragseries retrieved")

    if writefragseries:
        if not type(fragseries) == pd.core.series.Series:
            return "Supply fragment series to write to file"
        writepickle(fragseries, writefragseries)
        print("fragseries writted to file: " + writefragseries)

    # Step 4: Prepare cleaned pandas dataframe for fragment series attachment (output: fdbm)

    if (
        not type(fdbm) == pd.core.frame.DataFrame
        and not type(fdb) == pd.core.frame.DataFrame
    ):
        if not type(sdbdc) == dd.core.DataFrame:
            return "Please include a cleaned dask substance dataframe to which fragment information can be attached"
        fdbm = sdbdc.compute()
    #         fdbm.reset_index(inplace=True)
    #         fdbm.drop('index',axis=1,inplace=True)

    # Step 5: Attaching fragment series, generating an exploded fragment database (output: fdb)

    if type(fdbm) == pd.core.frame.DataFrame:
        if not type(fragseries) == pd.core.series.Series:
            return "Please include the fragment series that should be attached"
        fragdf = pd.DataFrame(fragseries)
        fdb = joindf(fragdf, fdbm, explodeDB=fragdf.columns[0])
        print("Unindexed fragment database completed")

    # Step 6: Adding additional columns, formatting, indexing and writing fragment database to file

    if fragsmiles:
        if not type(fdb) == pd.core.frame.DataFrame:
            return "Please supply fragment database to analyze"
        if not type(dfdb) == dd.core.DataFrame:
            if fdb.index.name or fdb.index.names:
                fdb.reset_index(inplace=True)
            dfdb = dd.from_pandas(fdb, npartitions=181)
            dfdb = client.persist(dfdb)
        fragseries = dfdb.map_partitions(
            getsmilespartition, meta=("FragmentSmiles", "O")
        ).compute()
        if writefragseries:
            if not type(fragseries) == pd.core.series.Series:
                return "Supply fragment series to write to file"
            writepickle(fragseries, writefragseries)
            print("fragseries writted to file: " + writefragseries)
        smilesdf = pd.DataFrame(fragseries)
        fdb = joindf(smilesdf, fdb)

    if mixtures:
        if not type(fdb) == pd.core.frame.DataFrame:
            return "Please supply fragment database to analyze"
        if not type(dfdb) == dd.core.DataFrame:
            if fdb.index.name or fdb.index.names:
                fdb.reset_index(inplace=True)
            dfdb = dd.from_pandas(fdb, npartitions=181)
            dfdb = client.persist(dfdb)
        fragseries = dfdb.map_partitions(
            findMixturespartition, meta=(">1 Compound", "boolean")
        )
        if writefragseries:
            if not type(fragseries) == pd.core.series.Series:
                return "Supply fragment series to write to file"
            writepickle(fragseries, writefragseries)
            print("fragseries writted to file: " + writefragseries)
        mixturedf = pd.DataFrame(fragseries)
        fdb = joindf(mixturedf, fdb)

    if mixturefrags:
        # NOTE: If exploded, fragment database should be aggregated by Smiles/Smarts
        if (
            not type(fdb) == pd.core.frame.DataFrame
            or "> 1 Compound" not in fdb.columns
        ):
            return "Please supply fragment database to analyze, with mixture indication (specify mixtures=True)"
        if not natoms:
            return "Please specify size of fragments retrieved"
        if not type(dfdb) == dd.core.DataFrame:
            dfdb = dd.from_pandas(fdb, npartitions=16)
            dfdb = client.persist(dfdb)
        fragseries = dfdb.map_partitions(
            getMixturefragspartition, natoms=natoms, meta=("FragmentSmarts", "O")
        )
        if writefragseries:
            if not type(fragseries) == pd.core.series.Series:
                return "Supply fragment series to write to file"
            writepickle(fragseries, writefragseries)
            print("fragseries writted to file: " + writefragseries)
        mixturesmarts = pd.DataFrame(fragseries)
        fdb = joindf(mixturesmarts, fdb)
    #%% Indexing dataframe

    if index:
        if not type(fdb) == pd.core.frame.DataFrame:
            return "Please supply fragment database to index"
        if fdb.index.name == index or fdb.index.names == index:
            print("Database is already indexed.")
        elif fdb.index.name or fdb.index.names:
            fdb.reset_index(inplace=True)
            fdb.set_index([index], inplace=True)
        else:
            fdb.set_index([index], inplace=True)
        print("Pandas fragment database indexed")
    if writefdb:
        writepickle(fdb, writefdb)
        print("Pandas fragment database written to file: " + writefdb)

    return fdb
