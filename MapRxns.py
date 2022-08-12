# %load ./MapRxns.py

from collections import (
    Counter,
    OrderedDict,
)
from math import ceil
import pandas as pd  # Preserve order of keys relative to reaction from left to right
from rdkit import Chem
from rdkit.Chem import rdChemReactions  # Reaction processing
import copy
from AnalgRxns import userefrxns
from rxnmapper import RXNMapper
from pydantic.typing import List
from BalanceRxns import balancerxn
from MainFunctions import (
    initray,
    molfromsmiles,
)  # Importing RXNMapper for unsupervised atom mapping
from FindFunctionalGroups import identify_functional_groups as IFG
import modin.pandas as mpd

#%% Reaction Mapping


def maprxn(rxns: List[str]):
    """
    For a given list of reactions, rxns, returns mapped reactions with confidence scores.
    Uses IBM transformer model.

    Parameters
    ----------
    rxns : List[str]
        List of reaction SMILES (no reactant/reagent split)

    Returns
    -------
    Output : list
        Mapped reactions with confidence scores

           mapped_rxn: str
               Mapped reaction SMARTS

           confidence: str
               Model confidence in the mapping rxn

    ['Error']: list
        If code doesn't run or mapper doesn't work
    """

    rxn_mapper = RXNMapper()
    try:
        return rxn_mapper.get_attention_guided_atom_maps(rxns)
    except Exception:
        return ["Error"]


def maprxns(row):
    """
    Applies maprxn to each row of a dataframe

    """
    balrxnsmiles = ""
    balrxnsmiles = row["balrxnsmiles"]
    mappedrxn = maprxn([balrxnsmiles])[0]
    if mappedrxn == "Error":
        return "Error", "Error"
    else:
        mapped_rxn = mappedrxn.get("mapped_rxn")
        conf = mappedrxn.get("confidence")
    return mapped_rxn, conf


def map_rxns(
    analoguerxnsbalfilt, refmappedrxns=None, ncpus=16, restart=True, reaxys_update=True
):  # Done
    """
    Applies maprxn to a given dataframe

    """
    #     breakpoint()
    if not analoguerxnsbalfilt.index.name and not analoguerxnsbalfilt.index.names:
        idxreset = True
    else:
        idxreset = False
    idxcol = []
    if reaxys_update:
        idxcol = ["ReactionID", "Instance"]
    else:
        idxcol = ["ReactionID"]
    if refmappedrxns is not None:
        analoguerxnsbalfilt, commondf = userefrxns(
            analoguerxnsbalfilt, idxcol=idxcol, refanaloguerxns=refmappedrxns
        )
        idxreset = False
    if not analoguerxnsbalfilt.empty:
        if ncpus > 1:
            if restart:
                initray(num_cpus=ncpus)
            if not idxreset:
                analoguerxnsbalfilt.reset_index(inplace=True)
                idxreset = True
            analoguerxnsbalfiltdis = mpd.DataFrame(analoguerxnsbalfilt)
        else:
            analoguerxnsbalfiltdis = analoguerxnsbalfilt
        mappedrxns = analoguerxnsbalfiltdis.apply(maprxns, axis=1, result_type="reduce")
        mappedrxns = pd.DataFrame(
            data=mappedrxns.values, index=mappedrxns.index, columns=["mappedrxns"]
        )
        mappedrxns[["mapped_rxn", "confidence"]] = pd.DataFrame(
            mappedrxns["mappedrxns"].tolist(), index=mappedrxns.index
        )
        mappedrxns.drop(columns=["mappedrxns"], inplace=True)
        analoguerxnsmapped = copy.deepcopy(analoguerxnsbalfilt)
        analoguerxnsmapped[["mapped_rxn", "confidence"]] = mappedrxns
        if idxreset:
            analoguerxnsmapped.set_index(idxcol, inplace=True)
        if refmappedrxns is not None and not commondf.empty:
            analoguerxnsmapped = pd.concat([analoguerxnsmapped, commondf])
    else:
        analoguerxnsmapped = commondf
    return analoguerxnsmapped


def checkrxns(
    analoguerxnsmappedfilt,
    refparsedrxns=None,
    ncpus=16,
    updateall=True,
    removeunmapped=True,
    restart=True,
    reaxys_update=True,
):  # Done
    if not analoguerxnsmappedfilt.index.name and not analoguerxnsmappedfilt.index.names:
        idxreset = True
    else:
        idxreset = False
    idxcol = []
    if reaxys_update:
        idxcol = ["ReactionID", "Instance"]
    else:
        idxcol = ["ReactionID"]
    if refparsedrxns is not None:
        analoguerxnsmappedfilt, commondf = userefrxns(
            analoguerxnsmappedfilt, idxcol=idxcol, refanaloguerxns=refparsedrxns
        )
        idxreset = False
    if not analoguerxnsmappedfilt.empty:
        if ncpus > 1:
            if restart:
                initray(num_cpus=ncpus)
            if not idxreset:
                analoguerxnsmappedfilt.reset_index(inplace=True)
                idxreset = True
            analoguerxnsmappedfiltdis = mpd.DataFrame(analoguerxnsmappedfilt)
        else:
            analoguerxnsmappedfiltdis = analoguerxnsmappedfilt
        compdupdate = analoguerxnsmappedfiltdis.apply(
            checkrxnrow,
            updateall=updateall,
            removeunmapped=removeunmapped,
            axis=1,
            result_type="reduce",
        )
        compdupdate = pd.Series(
            data=compdupdate.values, index=compdupdate.index
        )  # Optional convert modin back to pandas
        compdupdatedf = pd.DataFrame(
            data=compdupdate.tolist(),
            index=compdupdate.index,
            columns=["LHSdata", "RHSdata", "msg1"],
        )
        analoguerxnsparsed = copy.deepcopy(analoguerxnsmappedfilt)
        analoguerxnsparsed[["LHSdata", "RHSdata", "msg1"]] = compdupdatedf
        if idxreset:
            analoguerxnsparsed.set_index(idxcol, inplace=True)
        if refparsedrxns is not None and not commondf.empty:
            analoguerxnsparsed = pd.concat([analoguerxnsparsed, commondf])
    else:
        analoguerxnsparsed = commondf
    return analoguerxnsparsed


def checkrxnrow(row, updateall=True, removeunmapped=True):
    #     breakpoint()
    mappedrxn = row["mapped_rxn"]
    Rdata = row["LHSdata"]
    Pdata = row["RHSdata"]
    msg = row["msg"]
    if "with species" in msg:
        mandrcts = set(Rdata.keys()) - set(
            [
                int(addedspec)
                for addedspec in msg.rsplit("with species: ", 1)[1]
                .split(" with help product(s): ")[0]
                .split(",")
            ]
        )
    else:
        mandrcts = set(Rdata.keys())
    if "With hydrogen carriers" in msg:
        hcarriers = [
            int(hcarrier)
            for hcarrier in msg.split("With hydrogen carriers: ")[1]
            .split(", ")[0]
            .split(",")
        ]
    else:
        hcarriers = []
    if "Mandatory" in msg:
        mandrcts = mandrcts.union(
            {
                int(mandrct)
                for mandrct in msg.split("Mandatory species unmapped from LHS: ")[1]
                .split(", ")[0]
                .split(",")
            }
        )
    res = checkrxn(
        mappedrxn,
        Rdata=Rdata,
        Pdata=Pdata,
        updateall=updateall,
        removeunmapped=removeunmapped,
        mandrcts=mandrcts,
        hcarriers=hcarriers,
    )

    return res


def checkrxn(
    mappedrxn,
    Rdata={},
    Pdata={},
    ordered=True,
    updateall=True,
    removeunmapped=True,
    mandrcts=[],
    mandprods=[],
    hcarriers=[],
):  # Assume same rxn smiles stored next to each other
    """
    Checks reaction, updating mapped and clean molecules, removing unmapped species

    """

    #     breakpoint()
    rdrxn = rdChemReactions.ReactionFromSmarts(mappedrxn, useSmiles=True)
    cleanrxn = copy.copy(rdrxn)
    rdChemReactions.RemoveMappingNumbersFromReactions(cleanrxn)
    if ordered:
        LHSdata = OrderedDict({})
        RHSdata = OrderedDict({})
    else:
        LHSdata = {}
        RHSdata = {}
    msgr = []
    msgp = []
    # Updating LHSdata
    if Rdata:
        rsmilesfreq = gensmilesfreq(Rdata)
        smilesmismatch = []
        rmixtures = {}
        mismatch = False
        for ID, rct in enumerate(cleanrxn.GetReactants()):
            mappedmol = rdrxn.GetReactants()[ID]
            formula = Chem.rdMolDescriptors.CalcMolFormula(rct)
            if removeunmapped:
                if (
                    any(
                        [
                            atom.HasProp("molAtomMapNumber")
                            for atom in mappedmol.GetAtoms()
                        ]
                    )
                    or formula == "H2"
                    or hcarriers
                ):  # Confirmed, mapped reactant
                    LHSdata, rmixtures, msg_, ID0 = updatespecdict(
                        Rdata,
                        rsmilesfreq,
                        rct,
                        mappedmol,
                        updateddict=LHSdata,
                        mixtures=rmixtures,
                        updateall=updateall,
                        hcarriers=hcarriers,
                    )
                    if msg_ != "Valid":
                        #                         if ID0 not in smilesmismatch:
                        #                             smilesmismatch+=[ID0]
                        mismatch = True
            else:
                LHSdata, rmixtures, msg_, ID0 = updatespecdict(
                    Rdata,
                    rsmilesfreq,
                    rct,
                    mappedmol,
                    updateddict=LHSdata,
                    mixtures=rmixtures,
                    updateall=updateall,
                    hcarriers=hcarriers,
                )
                if msg_ != "Valid":
                    #                     if ID0 not in smilesmismatch:
                    #                         smilesmismatch+=[ID0]
                    mismatch = True
        if mismatch:
            smilesmismatch = [
                ID0 for ID0 in Rdata if ID0 not in LHSdata and ID0 not in rmixtures
            ]
            if smilesmismatch:
                msgr += [
                    "Smiles discrepancy for LHS species: "
                    + ", ".join([str(ID) for ID in smilesmismatch])
                ]
            else:
                msgr += ["Smiles discrepancy for LHS species"]
        #         breakpoint()
        if rmixtures:
            msgr += [
                "Mixture detected for LHS species: "
                + ",".join([str(ID) for ID in rmixtures])
            ]
            for ID0 in rmixtures:
                localcounts = [
                    int(Counter(rsmilesfreq[mixsmiles])[ID0])
                    for mixsmiles in rmixtures[ID0]
                ]
                numinsts = [
                    len(rmixtures[ID0][mixsmiles]) for mixsmiles in rmixtures[ID0]
                ]
                lb = [0 for j in range(len(numinsts))]
                ub = [0 for j in range(len(numinsts))]
                div = [
                    int(ceil(numinst / localcount))
                    for numinst, localcount in zip(numinsts, localcounts)
                ]
                count = max(div)
                if updateall:
                    LHSdata[ID0].update(
                        {"mappedsmiles": [], "cleanmol": [], "unmappedmix": []}
                    )
                    for i in range(count):
                        mappedsmiles = []
                        unmappedsmiles = []
                        cleanmol = []
                        for j, mixsmiles in enumerate(rmixtures[ID0]):
                            ub_ = min(numinsts[j], localcounts[j])
                            ub[j] = ub_ + lb[j]
                            mappedlist = rmixtures[ID0][mixsmiles][lb[j] : ub[j]]
                            mappedsmiles += [comb[0] for comb in mappedlist]
                            cleanmol += [comb[1] for comb in mappedlist]
                        mappedsmiles = tuple(mappedsmiles)
                        cleanmol = tuple(cleanmol)
                        unmappedsmiles = tuple(
                            [
                                mixsmiles
                                for mixsmiles in rsmilesfreq
                                for k in range(
                                    int(Counter(rsmilesfreq[mixsmiles])[ID0])
                                )
                                if mixsmiles in Rdata[ID0]["smiles"]
                                if mixsmiles not in rmixtures[ID0]
                            ]
                        )
                        LHSdata[ID0]["mappedsmiles"].extend([mappedsmiles])
                        LHSdata[ID0]["cleanmol"].extend([cleanmol])
                        LHSdata[ID0]["unmappedmix"].extend([unmappedsmiles])
                        lb = copy.deepcopy(ub)
                LHSdata[ID0]["count"] = count
        if removeunmapped:
            rem = Counter()
            Rspecies = Counter(
                [ID0 for ID0 in Rdata for _ in range(Rdata[ID0]["count"])]
            )
            LHSspecies = Counter(
                [ID0 for ID0 in LHSdata for _ in range(LHSdata[ID0]["count"])]
            )
            rem.update(Rspecies)
            rem.subtract(LHSspecies)
            #             removedmandrcts=[ID0 for ID0 in rem if ID0 in mandrcts if ID0 not in LHSspecies]
            #             removedrcts=[ID0 for ID0 in rem if ID0 not in removedmandrcts if ID0 not in smilesmismatch if ID0 in Rspecies if ID0 not in LHSspecies]
            removedrcts = [
                ID0
                for ID0 in rem
                if ID0 in Rspecies
                if ID0 not in LHSspecies
                if ID0 not in smilesmismatch
            ]
            removedmandrcts = [ID0 for ID0 in removedrcts if ID0 in mandrcts]
            # New
            removedmandrcts = list(
                set(removedmandrcts).union(
                    {mandrct for mandrct in mandrcts if mandrct not in LHSspecies}
                )
            )
            # New
            runmapped = [
                ID0
                for ID0 in rem
                if rem[ID0] > 0
                if ID0 not in removedrcts
                if ID0 not in smilesmismatch
            ]  #!=0
            if removedmandrcts:
                msgr += [
                    "Mandatory species unmapped from LHS: "
                    + ", ".join([str(ID) for ID in removedmandrcts])
                ]
                removedrcts = [ID0 for ID0 in removedrcts if ID0 not in removedmandrcts]
            if removedrcts:
                msgr += [
                    "Unmapped species from LHS: "
                    + ", ".join([str(ID) for ID in removedrcts])
                ]
            if runmapped:
                msgr += [
                    "Unmapped species instances from LHS: "
                    + ", ".join([str(ID) for ID in runmapped])
                ]
    if Pdata:
        #         breakpoint()
        psmilesfreq = gensmilesfreq(Pdata)
        smilesmismatch = []
        pmixtures = {}
        mismatch = False
        for ID, prod in enumerate(cleanrxn.GetProducts()):
            mappedmol = rdrxn.GetProducts()[ID]
            #             formula=rdkit.Chem.rdMolDescriptors.CalcMolFormula(prod)
            if removeunmapped:
                if any(
                    [atom.HasProp("molAtomMapNumber") for atom in mappedmol.GetAtoms()]
                ):  # Confirmed, mapped reactant
                    RHSdata, pmixtures, msg_, ID0 = updatespecdict(
                        Pdata,
                        psmilesfreq,
                        prod,
                        mappedmol,
                        updateddict=RHSdata,
                        mixtures=pmixtures,
                        updateall=updateall,
                    )
                    if msg_ != "Valid":
                        #                         if ID0 not in smilesmismatch:
                        #                             smilesmismatch+=[ID0]
                        mismatch = True
            else:
                RHSdata, pmixtures, msg_, ID0 = updatespecdict(
                    Pdata,
                    psmilesfreq,
                    prod,
                    mappedmol,
                    updateddict=RHSdata,
                    mixtures=pmixtures,
                    updateall=updateall,
                )
                if msg_ != "Valid":
                    #                     if ID0 not in smilesmismatch:
                    #                         smilesmismatch+=[ID0]
                    mismatch = True
        if mismatch:
            smilesmismatch = [
                ID0 for ID0 in Pdata if ID0 not in RHSdata and ID0 not in pmixtures
            ]
            if smilesmismatch:
                msgp += [
                    "Smiles discrepancy for RHS species: "
                    + ", ".join([str(ID) for ID in smilesmismatch])
                ]
            else:
                msgp += ["Smiles discrepancy for RHS species"]
        if pmixtures:
            msgp += [
                "Mixture detected for RHS species: "
                + ",".join([str(ID) for ID in pmixtures])
            ]
            #             breakpoint()
            for ID0 in pmixtures:
                localcounts = [
                    int(Counter(psmilesfreq[mixsmiles])[ID0])
                    for mixsmiles in pmixtures[ID0]
                ]
                numinsts = [
                    len(pmixtures[ID0][mixsmiles]) for mixsmiles in pmixtures[ID0]
                ]
                lb = [0 for j in range(len(numinsts))]
                ub = [0 for j in range(len(numinsts))]
                div = [
                    int(ceil(numinst / localcount))
                    for numinst, localcount in zip(numinsts, localcounts)
                ]
                count = max(div)
                if updateall:
                    RHSdata[ID0].update(
                        {"mappedsmiles": [], "cleanmol": [], "unmappedmix": []}
                    )
                    for i in range(count):
                        mappedsmiles = []
                        unmappedsmiles = []
                        cleanmol = []
                        for j, mixsmiles in enumerate(pmixtures[ID0]):
                            ub_ = min(numinsts[j], localcounts[j])
                            ub[j] = ub_ + lb[j]
                            mappedlist = pmixtures[ID0][mixsmiles][lb[j] : ub[j]]
                            mappedsmiles += [comb[0] for comb in mappedlist]
                            cleanmol += [comb[1] for comb in mappedlist]
                        mappedsmiles = tuple(mappedsmiles)
                        cleanmol = tuple(cleanmol)
                        unmappedsmiles = tuple(
                            [
                                mixsmiles
                                for mixsmiles in psmilesfreq
                                for k in range(
                                    int(Counter(psmilesfreq[mixsmiles])[ID0])
                                )
                                if mixsmiles in Pdata[ID0]["smiles"]
                                if mixsmiles not in pmixtures[ID0]
                            ]
                        )
                        RHSdata[ID0]["mappedsmiles"].extend([mappedsmiles])
                        RHSdata[ID0]["cleanmol"].extend([cleanmol])
                        RHSdata[ID0]["unmappedmix"].extend([unmappedsmiles])
                        lb = copy.deepcopy(ub)
                RHSdata[ID0]["count"] = count
        if removeunmapped:
            rem = Counter()
            Pspecies = Counter(
                [ID0 for ID0 in Pdata for _ in range(Pdata[ID0]["count"])]
            )
            RHSspecies = Counter(
                [ID0 for ID0 in RHSdata for _ in range(RHSdata[ID0]["count"])]
            )
            rem.update(Pspecies)
            rem.subtract(RHSspecies)
            #             removedmandprod=[ID0 for ID0 in rem if ID0 in mandprods if ID0 not in RHSspecies]
            #             removedprods=[ID0 for ID0 in rem if ID0 not in removedmandprod if ID0 not in smilesmismatch if ID0 in Pspecies if ID0 not in RHSspecies]
            removedprods = [
                ID0
                for ID0 in rem
                if ID0 in Pspecies
                if ID0 not in RHSspecies
                if ID0 not in smilesmismatch
            ]
            removedmandprods = [ID0 for ID0 in removedprods if ID0 in mandprods]
            punmapped = [
                ID0
                for ID0 in rem
                if rem[ID0] > 0
                if ID0 not in removedprods
                if ID0 not in smilesmismatch
            ]  #!=0
            if removedmandprods:
                msgp += [
                    "Mandatory species unmapped from RHS: "
                    + ", ".join([str(ID) for ID in removedmandprods])
                ]
                removedprods = [
                    ID0 for ID0 in removedprods if ID0 not in removedmandprods
                ]
            if removedprods:
                msgp += [
                    "Unmapped species from RHS: "
                    + ", ".join([str(ID) for ID in removedprods])
                ]
            if punmapped:
                msgp += [
                    "Unmapped species instances from RHS: "
                    + ", ".join([str(ID) for ID in punmapped])
                ]
    msg = msgr + msgp
    if not msg:
        msg = "Valid"
    else:
        msg = ", ".join(msg)
    return LHSdata, RHSdata, msg


def gensmilesfreq(specdict, validate=True):
    """
    Generates smile frequency dictionary (Sometimes species have the same SMILES with different IDs eg. mixture vs pure)

    """
    smilesfreq = {}
    for ID0 in specdict:
        specsmiles = specdict[ID0]["smiles"]
        specsmiles = specsmiles.split(".")
        for specsmile in specsmiles:
            if validate:
                specsmile = Chem.MolToSmiles(molfromsmiles(specsmile))
            if specsmile in smilesfreq:
                smilesfreq[specsmile].extend([ID0])
            else:
                smilesfreq.update({specsmile: [ID0]})
    return smilesfreq


def updatespecdict(
    refdict,
    smilesfreq,
    cleanmol,
    mappedmol,
    updateddict=OrderedDict({}),
    mixtures={},
    hcarriers=[],
    updateall=True,
):
    """
    Updates species dictionary based on given reactant and cleaned molecule from a reaction
    hcarriers (list of hydrogen containing species involved in reaction but not mapped)
    """
    foundmatch = False
    specsmiles = Chem.MolToSmiles(
        molfromsmiles(Chem.MolToSmiles(cleanmol))
    )  # Ensuring RDKit smiles
    #     breakpoint()
    if specsmiles not in smilesfreq:
        ID0 = ""
        msg = "Smiles discrepancy for species"
    else:
        idx = ""
        msg = "Valid"
        IDlist = smilesfreq[specsmiles]
        mixtures_ = ["." in refdict[ID0]["smiles"] for ID0 in IDlist]
        if len(IDlist) > 1:  # Try mixtures first
            pure = [
                ID0 for i, ID0 in enumerate(IDlist) if not mixtures_[i]
            ]  # Pure matches
            if any(mixtures):
                for i, ID0 in enumerate(IDlist):
                    if mixtures_[i] and ID0 in mixtures:
                        if specsmiles in mixtures[ID0]:
                            loccount = len(mixtures[ID0][specsmiles])
                        else:
                            loccount = 0
                        if any(
                            [
                                len(mixtures[ID0][specsmiles_]) > loccount
                                for specsmiles_ in mixtures[ID0]
                            ]
                        ):
                            idx = i
                            break
                if not idx and not pure:
                    for i, ID0 in enumerate(IDlist):
                        if mixtures_[i] and ID0 not in mixtures:
                            idx = i
                            break
            if not idx and pure:
                for i, ID0 in enumerate(IDlist):
                    if not mixtures_[i] and ID0 not in updateddict:
                        idx = i
                        break
            if not idx:
                idx = 0
        else:
            idx = 0
        ID0 = IDlist[idx]
        # New
        if (
            hcarriers
            and ID0 not in hcarriers
            and not any(
                [atom.HasProp("molAtomMapNumber") for atom in mappedmol.GetAtoms()]
            )
        ):
            return updateddict, mixtures, msg, ID0
        # New
        mixture = mixtures_[idx]
        if updateall:
            mappedsmiles = Chem.MolToSmiles(mappedmol)
            if mixture:
                if ID0 not in mixtures:
                    updateddict.update({ID0: copy.deepcopy(refdict[ID0])})
                    updateddict[ID0]["mixture"] = mixture
                    mixtures.update({ID0: {specsmiles: [(mappedsmiles, cleanmol)]}})
                elif specsmiles not in mixtures[ID0]:
                    mixtures[ID0].update({specsmiles: [(mappedsmiles, cleanmol)]})
                else:
                    mixtures[ID0][specsmiles].extend([(mappedsmiles, cleanmol)])
            else:
                if ID0 not in updateddict:
                    updateddict.update({ID0: copy.deepcopy(refdict[ID0])})
                    updateddict[ID0]["mixture"] = mixture
                    updateddict[ID0]["count"] = 1
                    updateddict[ID0].update(
                        {"mappedsmiles": [mappedsmiles], "cleanmol": [cleanmol]}
                    )
                else:
                    updateddict[ID0]["count"] += 1
                    updateddict[ID0]["mappedsmiles"].extend([mappedsmiles])
                    updateddict[ID0]["cleanmol"].extend([cleanmol])
        else:
            if mixture:
                if ID0 not in mixtures:
                    updateddict.update({ID0: copy.deepcopy(refdict[ID0])})
                    updateddict[ID0]["mixture"] = mixture
                    mixtures.update({ID0: {specsmiles: [()]}})
                elif specsmiles not in mixtures[ID0]:
                    mixtures[ID0].update({specsmiles: [()]})
                else:
                    mixtures[ID0][specsmiles].extend([()])
            else:
                if ID0 not in updateddict:
                    updateddict.update({ID0: copy.deepcopy(refdict[ID0])})
                    updateddict[ID0]["mixture"] = mixture
                    updateddict[ID0]["count"] = 1
                else:
                    updateddict[ID0]["count"] += 1
    return updateddict, mixtures, msg, ID0


def updaterxns(
    analoguerxnsparsed, hc_prod={}, analoguerxns=None, ncpus=16, restart=True
):
    """
    Updates reactions if there are unmapped species and balances if there are changes (optional)

    """
    if analoguerxns is not None:
        analoguerxnsparsed = updatecolumns(
            analoguerxns, analoguerxnsparsed, cols=["Rdata", "Rgtdata", "Solvdata"]
        )
    if ncpus > 1:
        if restart:
            initray(num_cpus=ncpus)
        analoguerxnsparseddis = mpd.DataFrame(analoguerxnsparsed)
    else:
        analoguerxnsparseddis = analoguerxnsparsed
    updatedrxn = analoguerxnsparseddis.apply(
        updaterxns_, hc_prod=hc_prod, axis=1, result_type="reduce"
    )
    updatedrxn = pd.DataFrame(
        data=updatedrxn.values, index=updatedrxn.index, columns=["rxncomb"]
    )
    analoguerxnsparsed[
        [
            "mapped_rxn",
            "confidence",
            "balrxnsmiles",
            "msg",
            "LHS",
            "RHS",
            "hcrct",
            "hcprod",
            "LHSdata",
            "RHSdata",
            "msg1",
        ]
    ] = pd.DataFrame(updatedrxn["rxncomb"].tolist(), index=updatedrxn.index)
    return analoguerxnsparsed


def updaterxns_(row, hc_prod={}):
    """
    Updates reactions if there are unmapped species and balances if there are changes (optional). Assumes both
    balancerxn, maprxn and checkrxn have all been called already

    """
    #     breakpoint()
    msg1 = copy.deepcopy(row["msg1"])
    msg = copy.deepcopy(row["msg"])  # Balanced output message
    LHSdata = copy.deepcopy(row["LHSdata"])
    RHSdata = copy.deepcopy(row["RHSdata"])
    hcprod = copy.deepcopy(row["hcprod"])
    hcrct = copy.deepcopy(row["hcrct"])
    if "Rgtdata" in row.keys():
        Rgtdata = row["Rgtdata"]
    else:
        Rgtdata = {}
    if "Solvdata" in row.keys():
        Solvdata = row["Solvdata"]
    else:
        Solvdata = {}
    if "Rdata" in row.keys():
        mandrcts = row["Rdata"]
    else:
        mandrcts = LHSdata
    addedspecies = list(set(LHSdata.keys()) - set(mandrcts.keys()))
    #     else:
    #         if 'with species' in msg:
    #             mandrcts=set(LHSdata.keys())-{int(addedspec) for addedspec in msg.split('with species: ')[1].split(' with help product(s): ')[0].split(',')}
    #         else:
    #             mandrcts=set(LHS.keys())

    #         if 'Mandatory' in msg1:
    #             mandrcts=mandrcts.union({int(mandrct) for mandrct in msg.split('Mandatory species unmapped from LHS: ')[1].split(', ')[0].split(',')})
    #     breakpoint()
    storemsg = ""
    i = 0
    while "Unmapped" in msg1 or i == 0:  # Unmapped species exist not reflected
        if "from RHS" in msg1:  # Mandatory products unmapped
            storemsg = msg1
        if "Smiles discrepancy" in msg1:
            break
        if hcprod is not None:
            hcprod = [hcprod_ for hcprod_ in hcprod if hcprod_ in RHSdata]
        hcrct = [hcrct_ for hcrct_ in hcrct if hcrct_ in LHSdata]
        balrxnsmiles, msg, LHS, RHS, hcrct, hcprod, LHSdata, RHSdata = balancerxn(
            LHSdata,
            RHSdata,
            first=False,
            Rgtdata=Rgtdata,
            Solvdata=Solvdata,
            addedspecies=addedspecies,
            hc_prod=hc_prod,
            coefflim=6,
            mandrcts=mandrcts,
            usemapper=False,
            ignoreH=False,
        )
        #         if 'Hydrogen carriers' in msg or not hc_prod: #No point balancing again as hydrogen deficit always present
        #             balrxnsmiles,_,LHSids,RHSids,_,_,_,_=update_rxn(LHSdata,RHSdata,hc_prod=hc_prod,hcprod=hcprod,hcrct=hcrct,msg=msg)
        mappedrxn = maprxn([balrxnsmiles])[0]
        if mappedrxn == "Error":
            mapped_rxn = "Error"
            conf = "Error"
            msg1 = "Mapping error"
            break
        else:
            mapped_rxn = mappedrxn.get("mapped_rxn")
            conf = mappedrxn.get("confidence")
            if "With hydrogen carriers" in msg:
                hcarriers = [
                    int(hcarrier)
                    for hcarrier in msg.split("With hydrogen carriers: ")[1]
                    .split(", ")[0]
                    .split(",")
                ]
            else:
                hcarriers = []
            LHSdata, RHSdata, msg1 = checkrxn(
                mapped_rxn,
                Rdata=LHSdata,
                Pdata=RHSdata,
                updateall=True,
                removeunmapped=True,
                mandrcts=mandrcts,
                hcarriers=hcarriers,
            )
        if storemsg:
            if msg1 == "Valid":
                msg1 = storemsg
            else:
                msg1 = storemsg + ", " + msg1
            break
        i += 1
    return (
        mapped_rxn,
        conf,
        balrxnsmiles,
        msg,
        LHS,
        RHS,
        hcrct,
        hcprod,
        LHSdata,
        RHSdata,
        msg1,
    )


def updatecolumns(parent, child, cols=[], config=[]):
    if type(parent) == str:
        parent = pd.read_pickle(parent)
    if type(child) == str:
        child = pd.read_pickle(child)
    child[cols] = copy.deepcopy(parent[cols])
    if config:
        child = child[config]
    return child


def assignfrags(analoguerxnsparsedfilt, fragdict, strict=False, ncpus=16, restart=True):
    """
    Assigns fragments to analogue species, and corresponding atom indices/mapping numbers


    """
    if ncpus > 1:
        if restart:
            initray(num_cpus=ncpus)
        analoguerxnsparsedfiltdis = mpd.DataFrame(analoguerxnsparsedfilt)
    else:
        analoguerxnsparsedfiltdis = analoguerxnsparsedfilt
    compdassigned = analoguerxnsparsedfiltdis.apply(
        assignfragsrow, fragdict=fragdict, strict=strict, axis=1, result_type="reduce"
    )
    compdassigned = pd.Series(
        data=compdassigned.values, index=compdassigned.index
    )  # Optional convert modin back to pandas
    compdassigneddf = pd.DataFrame(
        data=compdassigned.tolist(),
        index=compdassigned.index,
        columns=["LHSdata", "nofg", "msg2"],
    )
    analoguerxnsassigned = copy.deepcopy(analoguerxnsparsedfilt)
    analoguerxnsassigned[["LHSdata", "nofg", "msg2"]] = compdassigneddf
    return analoguerxnsassigned


def assignfragsrow(row, fragdict, strict=False):
    """
    Assigns fragments to analogue species, and corresponding atom indices/mapping numbers

    """
    LHSdata = row.LHSdata
    return assignfrags_(LHSdata, fragdict, strict=strict)


def assignfrags_(LHSdata, fragdict, strict=False):
    """
    Assigns fragments to analogue species, and corresponding atom indices/mapping numbers

    """
    if strict:
        checkresults = True  # Not recommended now as isotopic, hydrocarbon and user-defined fragments may be rejected
    else:
        checkresults = False
    fragdata = copy.deepcopy(LHSdata)
    msg = []
    nofg = set()  # Species with no functional groups
    nonanalogue = set()  # Species that don't have correct/desired functional group
    for rctid in fragdata:
        fragloc = {}
        fragdata[rctid]["querycompds"] = {}
        commonfrags = [
            frag for frag in fragdict if rctid in fragdict[frag]["analoguepool"]
        ]
        if not commonfrags:
            nonanalogue.add(rctid)
            continue
        fragdata[rctid]["querycompds"].update(
            {frag: fragdict[frag]["Query"] for frag in commonfrags}
        )
        for idx, cleanmol0 in enumerate(fragdata[rctid]["cleanmol"]):
            if type(cleanmol0) == tuple:  # Mixture detected
                for idx2, cleanmol in enumerate(cleanmol0):
                    cleanmol = Chem.AddHs(cleanmol)
                    for frag in commonfrags:
                        fragloc, nofg = update_matches(
                            cleanmol,
                            frag,
                            fragloc=fragloc,
                            nofg=nofg,
                            checkresults=checkresults,
                            idx=(idx, idx2),
                            rctid=rctid,
                        )
                    if not fragloc[(idx, idx2)]:
                        nonanalogue.add(
                            rctid
                        )  # This will remove all mixtures (some parts of the mixture may not react)
            else:
                cleanmol = Chem.AddHs(cleanmol0)  # Add hydrogens
                for frag in commonfrags:
                    fragloc, nofg = update_matches(
                        cleanmol,
                        frag,
                        fragloc=fragloc,
                        nofg=nofg,
                        checkresults=checkresults,
                        idx=idx,
                        rctid=rctid,
                    )
                if not fragloc[idx]:
                    nonanalogue.add(rctid)
        fragdata[rctid]["fragloc"] = fragloc
    if nonanalogue:
        msg = (
            "Species "
            + ", ".join([str(rctid) for rctid in nonanalogue])
            + " not analogue"
        )
    else:
        msg = "Valid"
    return fragdata, nofg, msg


#%% Substructure matching
def update_matches(
    mol, pattsmiles, checkresults=False, fragloc={}, nofg=set(), idx=0, rctid=0
):
    """
    Returns atom indices for substructure matches of a pattern in a molecule
    """
    #     breakpoint()
    patt = Chem.MolFromSmarts(pattsmiles)
    patt.UpdatePropertyCache(strict=False)
    corr_matches, funcgroupids, msg = get_matches(
        mol, patt, checkresults=checkresults
    )  # funcgroupids refers to active fragment, change checkresults to true if more strict
    if not corr_matches:
        fragloc.update({idx: {}})
    else:
        if msg != "Valid" or not funcgroupids:
            nofg.add(rctid)
        #         if not includeHs:
        #             corr_matches=[tuple((idx for idx in co if idx<len(Chem.RemoveHs(mol).GetAtoms()))) for co in corr_matches]
        if idx not in fragloc.keys():
            fragloc.update(
                {
                    idx: {
                        pattsmiles: {
                            "corrmatches": corr_matches,
                            "funcgroupids": funcgroupids,
                        }
                    }
                }
            )
        else:
            fragloc[idx].update(
                {
                    pattsmiles: {
                        "corrmatches": corr_matches,
                        "funcgroupids": funcgroupids,
                    }
                }
            )
    return fragloc, nofg


def get_matches(mol, patt, checkresults=False):
    """
    Returns atom indices for substructure matches of a pattern in a molecule

        Parameters
        ----------
        mol : RDKit mol
            Molecule to check pattern
        patt : RDKit mol
            Pattern fragment
        checkresults: bool
            Optional, True if strict match with the pattern is needed including hydrogens

        Returns
        -------
       corr_matches: Set
           Set of tuples containing atom indices for every pattern match, after verification
           that correct functional group is present. Returned only if checkresults is True.

       funcgroupids: List
           List of atom ids corresponding to functional groups in each correct pattern match

        matches: Set
            Set of tuples containing atom indices for every pattern match. Returned only if checkresults is False.
    """
    msg = ""
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return False, False, "No match"
    funcgroupmol = IFG(mol)  # Functional groups of RDKit reactant
    funcgrouppatt = IFG(patt)  # Functional groups of carrier fragment
    if not funcgroupmol:
        msg = "No functional group in parent"
        checkresults = False
    elif not funcgrouppatt:
        msg = "No functional group in pattern"
        checkresults = False
    if checkresults:  # Buggy
        #         breakpoint()
        funcids = (
            set()
        )  # Store functional groups that are of the same type as the carrier fragment
        for funcgroup in funcgrouppatt:
            matchtype = [
                molgroup
                for molgroup in funcgroupmol
                if molgroup.atoms == funcgroup.atoms
            ]  # change to .atoms if not working
            for molgroup in matchtype:
                #                 if not any([atoms_are_different(mol.GetAtomWithIdx(atomid),patt.GetAtomWithIdx(pattid),usesmarts=False)
                #                             for atomid,pattid in zip(molgroup.atomIds,funcgroup.atomIds)]): #BUGGY
                funcids.update({atomid for atomid in molgroup.atomIds})
        corr_matches = [match for match in matches if set(match).intersection(funcids)]
        funcgroupids = [set(match).intersection(funcids) for match in corr_matches]
        msg = "Valid"
    else:
        #         breakpoint()
        funcgroupids = []  # Added
        corr_matches = [match for match in matches]
        if not msg:
            funcids = {
                atomid for molgroup in funcgroupmol for atomid in molgroup.atomIds
            }
            funcgroupids = [set(match).intersection(funcids) for match in corr_matches]
            msg = "Valid"
    return corr_matches, funcgroupids, msg
