# %load ./ImpurityRanking.py

import re
from matplotlib import pyplot as plt
import pandas as pd
import modin.pandas as mpd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from collections import OrderedDict
import numpy as np
import cirpy
import copy
from math import ceil, floor
from AnalgCompds import getCarrierFrags0
from MainFunctions import getfragments, initray, molfromsmiles, openpickle
from MapRxns import get_matches


def updatecatalyst(impfinalfilt, catpool, ncpus=16, restart=True):
    if ncpus > 1:
        if restart:
            initray(num_cpus=ncpus)
        impfinalfiltdis = mpd.DataFrame(impfinalfilt)
    else:
        impfinalfiltdis = impfinalfilt
    updatedcatalyst = impfinalfiltdis.apply(
        updatecatalyst_, catpool=catpool, axis=1, result_type="reduce"
    )
    updatedcatalyst = pd.Series(
        data=updatedcatalyst.values, index=updatedcatalyst.index
    )
    impfinalfilt["CatalystID2"] = updatedcatalyst
    return impfinalfilt


def updatecatalyst_(row, catpool):
    #     breakpoint()
    rgt = row["ReagentID"]
    cat = row["CatalystID"]
    if rgt == "NaN" or rgt is None or not rgt:
        return cat
    else:
        combinedcat = list((set(rgt).intersection(catpool))) + cat
    return combinedcat


def relevance_score_morgan(
    impfinalfilt,
    fragdict,
    version="mean",
    includereagents=False,
    expand=1,
    morganradius=2,
    ncpus=16,
    restart=False,
):
    if ncpus > 1:
        if restart:
            initray(num_cpus=ncpus)
        impfinalfiltdis = mpd.DataFrame(impfinalfilt)
    else:
        impfinalfiltdis = impfinalfilt
    relevance_morgan = impfinalfiltdis.apply(
        relevance_score_morgan_,
        fragdict=fragdict,
        version=version,
        includereagents=includereagents,
        expand=expand,
        morganradius=morganradius,
        axis=1,
        result_type="reduce",
    )
    relevance_morgan = pd.Series(
        data=relevance_morgan.values, index=relevance_morgan.index
    )
    impfinalfilt["Relevance_morgan"] = relevance_morgan
    return impfinalfilt


def relevance_score_morgan_(
    row, fragdict, version="mean", includereagents=False, expand=1, morganradius=2
):  # Will use morgan/circular fingerprints and dice similarity
    #     breakpoint()
    LHSids = [
        compd for compd in row.LHSdata for _ in range(row.LHSdata[compd]["count"])
    ]
    querycompds = copy.deepcopy(row.querycompds)
    LHSdata = copy.deepcopy(row.LHSdata)
    Rgtdata = copy.deepcopy(row.Rgtdata)
    catids = copy.deepcopy(row.CatalystID2)
    pos = 0
    cycles = 0
    if version == "mean":
        globalsim = 0
    else:
        globalsim = 1
    #     breakpoint()
    while pos <= len(querycompds) - 1:
        cycles += 1
        query = querycompds[pos]
        analgcompd = LHSids[pos]
        analgsmiles = LHSdata[analgcompd]["smiles"]
        m1 = Chem.MolFromSmiles(query)
        m2 = Chem.MolFromSmiles(analgsmiles)  # 4-nitrophenol
        fp1 = AllChem.GetMorganFingerprint(m1, morganradius)
        fp2 = AllChem.GetMorganFingerprint(m2, morganradius)
        localsim = DataStructs.DiceSimilarity(fp1, fp2)
        if version == "mean":
            globalsim += localsim
        else:
            globalsim = globalsim * localsim
        pos += row.LHSdata[analgcompd]["count"]
    if includereagents and Rgtdata:
        rgtids = [
            rgtid
            for rgtid in list(Rgtdata.keys())
            if rgtid not in LHSids and rgtid not in catids
        ]
        for rgtid in rgtids:
            cycles += 1
            rgtsmiles = Rgtdata[rgtid]["smiles"]
            m1 = Chem.MolFromSmiles(rgtsmiles)
            fp1 = AllChem.GetMorganFingerprint(m1, morganradius)
            frags = set(
                getCarrierFrags0(m1, userinput="mol", expand=expand, resFormat="smiles")
            ).intersection(set(fragdict.keys()))
            if not frags:
                frags = [
                    frag for frag in fragdict if rgtid in fragdict[frag]["analoguepool"]
                ]
            if frags:
                querycompds = []
                for frag in frags:
                    querycompds += fragdict[frag]["Query"]
                localsim = []
                for querycompd in querycompds:
                    m2 = Chem.MolFromSmiles(querycompd)
                    fp2 = AllChem.GetMorganFingerprint(m2, 2)
                    localsim.append(DataStructs.DiceSimilarity(fp1, fp2))
                localsim = max(localsim)
            else:
                localsim = 0
            if version == "mean":
                globalsim += localsim
            else:
                globalsim = globalsim * localsim
    if version == "mean":
        return globalsim / cycles
    else:
        return globalsim


def standardize(impfinalfilt, reaxys_update=True, ncpus=16, restart=False):
    if ncpus > 1:
        if restart:
            initray(num_cpus=ncpus)
        impfinalfiltdis = mpd.DataFrame(impfinalfilt)
    else:
        impfinalfiltdis = impfinalfilt
    standardized = impfinalfilt.apply(standardizerxns, axis=1, result_type="reduce")
    standardized = pd.DataFrame(
        data=standardized.tolist(),
        index=standardized.index,
        columns=["querycompds", "impurities", "impurityrxn"],
    )
    impfinalfilt[["querycompds", "impurities", "impurityrxn"]] = standardized
    if reaxys_update:
        impfinalfilt = (
            impfinalfilt.reset_index()
            .drop_duplicates(subset=["Instance", "ReactionID", "impurityrxn"])
            .set_index(["ReactionID", "Instance"])
        )
    else:
        impfinalfilt = (
            impfinalfilt.reset_index()
            .drop_duplicates(subset=["ReactionID", "impurityrxn"])
            .set_index("ReactionID")
        )
    return impfinalfilt


def standardizerxns(row):
    querycompds = copy.deepcopy(row.querycompds)
    impurities = copy.deepcopy(row.impurities)
    querycompds = tuple(sorted(querycompds))
    impurities = tuple(sorted(impurities))
    impurityrxn = ">>".join(
        [getfragments(querycompds, smiles=True), getfragments(impurities, smiles=True)]
    )
    return querycompds, impurities, impurityrxn


# def translatenotes(row):
# #     breakpoint()
#     conditionnotes=copy.deepcopy(row.ConditionNotes)
#     if conditionnotes:
#         translator=Translator()
#         time.sleep(4)
# #         try:
#         translated=translator.translate(conditionnotes)
#         translated=translated.text
# #         except Exception:
# #             return conditionnotes
#         m = re.search(' <sub> (.+?) </ sub>', translated)
#         if m:
#             found = m.group(1)
#             re.sub(' <sub> (.+?) </ sub>','<sub>'+found+'</sub',translated)
#         return translated
#     else:
#         return conditionnotes

# def checkmissing_(row,analoguerxnsfinal):
#     rct=row.LHS
#     ID=row.name
#     refrct=analoguerxnsfinal.loc[ID].ReactantID
#     prod=row.RHS
#     refprod=analoguerxnsfinal.loc[ID].ProductID

#     if len(rct)<len(refrct) or len(prod)<len(refprod):
#         return False
#     else:
#         return True

# def checkmissing(impfinalfilt,analoguerxnssource):
#     if type(analoguerxnssource)==str:
#         analoguerxnsfinal=pd.read_pickle(analoguerxnssource)
#     elif type(analoguerxnssource)==pd.core.frame.DataFrame:
#         analoguerxnsfinal=analoguerxnssource
#     initray()
#     impfinalfiltdis=mpd.DataFrame(impfinalfilt)
#     mask=impfinalfiltdis.apply(checkmissing_,analoguerxnsfinal=analoguerxnsfinal,axis=1,result_type='reduce')
#     mask=pd.Series(data=mask.values,index=mask.index)
#     masklist=mask[mask.values==True].index
#     impfinalfilt=impfinalfilt[impfinalfilt.index.isin(masklist)]
#     return impfinalfilt


def conditionfilter(
    impfinalfilt, conditions=[]
):  # German entries are not translated/covered yet
    conditionlist = [
        "radiat",
        "pyrolysis",
        "enzym",
        "sonic",
        "electro",
        "distillation",
        "bacter",
        "gasif",
    ]
    conditionlist = [
        condition
        for condition in conditionlist
        if not any(
            [(condition in cond.lower()) or (condition in cond) for cond in conditions]
        )
    ]
    for condition in conditionlist:
        impfinalfilt = impfinalfilt.loc[
            ~impfinalfilt.ConditionNotes.str.contains(condition, case=False)
        ]
    return impfinalfilt


def catfilter_(row, catsmiles, catnames, substancedb, unresolvedids):
    #     breakpoint()
    catID = copy.deepcopy(row.CatalystID2)
    missingcat = copy.deepcopy(row.MissingCatalyst)
    namedict = copy.deepcopy(row.NameDict)
    if not catID and not missingcat:
        return True
    if catID:
        for catID_ in catID:
            matched = False
            for catsmiles_, catnames_ in zip(catsmiles, catnames):
                if (
                    catID_ in unresolvedids
                    or catsmiles_ is None
                    or catID_ not in substancedb
                ):  # No smiles
                    if namedict[catID_] is not None and catnames_ is not None:
                        if any(
                            [
                                re.match(
                                    catnames__, namedict[catID_], flags=re.IGNORECASE
                                )
                                for catnames__ in catnames_
                            ]
                        ):
                            matched = True
                            break
                        else:
                            continue
                    else:
                        return False
                else:
                    if substancedb.loc[catID_].Smiles == catsmiles_:
                        matched = True
                        break
                    else:
                        continue
            if not matched:
                return False

    if missingcat:
        for missingcat_ in missingcat:
            matched = False
            for catnames_ in catnames:
                if catnames_ is not None:
                    if any(
                        [
                            re.match(catnames__, missingcat_, flags=re.IGNORECASE)
                            for catnames__ in catnames_
                        ]
                    ):
                        matched = True
                        break
                    else:
                        continue
                else:
                    return False
            if not matched:
                return False
    return True


def catfilter(impfinalfilt, substancesource, unresolveddir, catalyst=[], useray=True):
    if not catalyst:
        impfinalfilt = impfinalfilt.loc[
            (~impfinalfilt.CatalystID2.astype(bool))
            & (~impfinalfilt.MissingCatalyst.astype(bool))
        ]
    else:
        if type(substancesource) == str:
            substancedb = pd.read_pickle(substancesource)
        elif type(substancesource) == pd.core.frame.DataFrame:
            substancedb = substancesource
        if type(unresolveddir) == str:
            unresolvedids = openpickle(unresolveddir)  # Unresolved IDs
        elif type(unresolveddir) == list:
            unresolvedids = unresolveddir
        catsmiles = []
        catnames = []
        for catref in catalyst:
            try:
                catmol = molfromsmiles(catref)
                catsmiles_ = Chem.MolToSmiles(catmol)
                smiles = True
            except Exception:
                smiles = False
            if smiles:
                catsmiles += [catsmiles_]
                catnames_ = cirpy.resolve(catsmiles_, "names")
                catnames_ = list(map(re.escape, catnames_))
                catnames += [catnames_]
            else:
                catnames_ = cirpy.resolve(catref, "names")
                catnames_ = list(map(re.escape, catnames_))
                catnames += [catnames_]
                catsmiles += [None]  # Can try to resolve but better not to avoid errors
        if useray:
            impfinalfiltdis = mpd.DataFrame(impfinalfilt)
            mask = impfinalfiltdis.apply(
                catfilter_,
                catsmiles=catsmiles,
                catnames=catnames,
                substancedb=substancedb,
                unresolvedids=unresolvedids,
                axis=1,
                result_type="reduce",
            )
            mask = pd.Series(data=mask.values, index=mask.index)
        else:
            mask = impfinalfilt.apply(
                catfilter_,
                catsmiles=catsmiles,
                catnames=catnames,
                substancedb=substancedb,
                unresolvedids=unresolvedids,
                axis=1,
                result_type="reduce",
            )
        masklist = mask[mask.values == True].index
        impfinalfilt = impfinalfilt[impfinalfilt.index.isin(masklist)]

    return impfinalfilt


def updatetemp(row):
    #     breakpoint()
    temp = copy.deepcopy(row.Temperature)
    updatedt = []
    if not temp:
        return updatedt
    if type(temp) != list:
        temp = [temp]
    for t in temp:
        t = str(t)
        t2 = [elem.strip() for elem in re.split(",|-", t)]
        t3 = []
        if not all(t2):  # Some negative temperatures exist
            negative = False
            for i, t2elem in enumerate(t2):
                if not t2elem:
                    negative = True
                    t3elem = "-"
                elif negative:
                    t3elem += t2elem
                    t3 += [t3elem]
                    negative = False
                else:
                    t3 += [t2elem]
        else:  # No negative temperatures
            t3 = t2
        t3 = [round(float(t3elem)) for t3elem in t3]
        maxtemp = max(t3)
        mintemp = min(t3)
        if maxtemp != mintemp:
            t3 = list(
                range(mintemp, maxtemp + 1, 1)
            )  # Note float values will be rounded to int values
        updatedt += t3
    return updatedt


# Fragment expansion and comparison


def popfragments(specieslist, expandfrag=[1, 2, 3]):
    fraginfo = {}
    for expand in expandfrag:
        fragloc = {}
        for queryspec in specieslist:
            m = Chem.AddHs(molfromsmiles(queryspec))
            if queryspec in fragloc:
                continue
            else:
                fragloc.update(OrderedDict({queryspec: OrderedDict()}))
            #     frag0dat=IFG(m)
            #     frag0=[grp.atoms for grp in frag0dat]
            #     frag0idx=[grp.atomIds for grp in frag0dat]
            #     fragdat.update({queryspec:{'frag0':frag0,'frag0idx':frag0idx}})
            fragdat = getCarrierFrags0(queryspec, expand=expand, resFormat="smiles")
            if type(fragdat) != list:
                fragdat = [fragdat]
            for i, fragsmiles in enumerate(fragdat):
                molfrag = ""
                if fragsmiles in fragloc[queryspec]:
                    continue
                pattmol = Chem.MolFromSmarts(fragsmiles)
                pattmol.UpdatePropertyCache(strict=False)
                corr_matches, funcids, _ = get_matches(
                    m, pattmol, checkresults=False
                )  # Chem.MolFromSmarts(Chem.MolToSmarts(Chem.MolFromSmiles(fragsmiles,sanitize=False)))
                if not corr_matches:
                    continue
                if expand > 1:  # Check if any other functional groups are included
                    #                     breakpoint()
                    fragdat1 = fraginfo["fraginfo1"][queryspec]
                    fraglist1 = list(fragdat1.keys())
                    if i > len(fraglist1) - 1:
                        i = len(fraglist1) - 1
                    funcid1 = fragdat1[fraglist1[i]]["funcgroupids"]
                    if funcid1 != funcids:
                        for j, func in enumerate(funcids):
                            if j > len(funcid1) - 1:
                                func1 = funcid1[len(funcid1) - 1]
                            else:
                                func1 = funcid1[j]
                            if func1 != func:
                                extra = func - func1
                                for funcid_1 in list(fragdat1.values()):
                                    if not extra:
                                        break
                                    for k, funcid__1 in enumerate(
                                        funcid_1["funcgroupids"]
                                    ):
                                        if not extra:
                                            break
                                        inters = funcid__1.intersection(extra)
                                        if inters:
                                            corr_matches[j] = tuple(
                                                set(corr_matches[j]).union(
                                                    set(funcid_1["corrmatches"][k])
                                                )
                                            )
                                            molfrag = Chem.MolFragmentToSmiles(
                                                m, corr_matches[j]
                                            )
                                            extra = extra - inters

                #                 corr_matches=[tuple((idx for idx in co if idx<len(molfromsmiles(queryspec).GetAtoms()))) for co in corr_matches]
                if not molfrag:
                    fragloc[queryspec].update(
                        OrderedDict(
                            {
                                fragsmiles: {
                                    "corrmatches": corr_matches,
                                    "funcgroupids": funcids,
                                }
                            }
                        )
                    )
                else:
                    fragloc[queryspec].update(
                        OrderedDict(
                            {
                                molfrag: {
                                    "corrmatches": corr_matches,
                                    "funcgroupids": funcids,
                                }
                            }
                        )
                    )

        fraginfo.update({"fraginfo" + str(expand): fragloc})
    return fraginfo


def popfragmentsrow(row, includereagents=True):
    balrxnsmiles = copy.deepcopy(row.balrxnsmiles)
    Rgtdata = copy.deepcopy(row.Rgtdata)
    catids = copy.deepcopy(row.CatalystID2)
    splitrxn = balrxnsmiles.split(">>")
    rcts = set(splitrxn[0].split("."))
    if includereagents and Rgtdata:
        rgtids = [rgtid for rgtid in list(Rgtdata.keys()) if rgtid not in catids]
        rgtsmiles = [Rgtdata[rgtid]["smiles"] for rgtid in rgtids]
        rcts = rcts.union(set(rgtsmiles))
    fraginfo = popfragments(rcts)
    return fraginfo


def fragcomparerow(row, queryfraginfoa, includefp=True):
    querycompds = copy.deepcopy(row.querycompds)
    LHSdata = copy.deepcopy(row.LHSdata)
    LHSsmiles = [LHSdata[rctid]["smiles"] for rctid in LHSdata]
    fragcompare = OrderedDict(
        {expand: OrderedDict(OrderedDict()) for expand in queryfraginfoa.keys()}
    )
    for i, expand in enumerate(queryfraginfoa.keys()):
        queryfraginfo = queryfraginfoa[expand]
        fraginfo = copy.deepcopy(row.Fraginfo[expand])
        if i == 0:  # Assuming fragment size of 1 used for analogue
            for j, query in enumerate(dict.fromkeys(querycompds).keys()):
                #                 breakpoint()
                if query in fragcompare[expand].keys():
                    continue
                if j > len(list(fraginfo.keys())) - 1:
                    j = len(list(fraginfo.keys())) - 1
                analg = list(fraginfo.keys())[j]
                commonfrag = set(queryfraginfo[query].keys()).intersection(
                    set(fraginfo[analg].keys())
                )
                if not commonfrag:
                    for k in range(len(fraginfo.keys())):
                        if k == j:
                            continue
                        analg = list(fraginfo.keys())[k]
                        commonfrag = set(queryfraginfo[query].keys()).intersection(
                            set(fraginfo[analg].keys())
                        )
                        if commonfrag:
                            break
                uniquefragsa = set(fraginfo[analg].keys()) - set(
                    queryfraginfo[query].keys()
                )
                uniquefragsq = set(queryfraginfo[query].keys()) - set(
                    fraginfo[analg].keys()
                )
                # Fingerprints
                if includefp:
                    m1 = molfromsmiles(query)
                    m2 = molfromsmiles(analg)
                    fp1m = AllChem.GetMorganFingerprint(m1, 2)
                    fp2m = AllChem.GetMorganFingerprint(m2, 2)
                    simm = DataStructs.DiceSimilarity(fp1m, fp2m)
                    fp1t = FingerprintMols.FingerprintMol(m1)
                    fp2t = FingerprintMols.FingerprintMol(m2)
                    simt = DataStructs.FingerprintSimilarity(fp1t, fp2t)
                if analg in LHSsmiles:
                    role = "rct"
                else:
                    role = "agent"
                for expand in queryfraginfoa.keys():
                    if query not in fragcompare[expand]:
                        fragcompare[expand].update(
                            OrderedDict(
                                {
                                    query: {
                                        analg: {
                                            "role": role,
                                            "morgansim": round(simm, 2),
                                            "topsim": round(simt, 2),
                                        }
                                    }
                                }
                            )
                        )
                    else:
                        fragcompare[expand][query].update(
                            {
                                analg: {
                                    "role": role,
                                    "morgansim": round(simm, 2),
                                    "topsim": round(simt, 2),
                                }
                            }
                        )
                expand = list(queryfraginfoa.keys())[i]
                fragcompare[expand][query][analg].update(
                    {
                        "commonfrags": list(commonfrag),
                        "uniquefragsa": list(uniquefragsa),
                        "uniquefragsq": list(uniquefragsq),
                    }
                )

        else:
            for query in fragcompare[expand]:
                for analg in fragcompare[expand][query]:
                    commonfrag = set(queryfraginfo[query].keys()).intersection(
                        set(fraginfo[analg].keys())
                    )
                    uniquefragsa = set(fraginfo[analg].keys()) - set(
                        queryfraginfo[query].keys()
                    )
                    uniquefragsq = set(queryfraginfo[query].keys()) - set(
                        fraginfo[analg].keys()
                    )
                    fragcompare[expand][query][analg].update(
                        {
                            "commonfrags": list(commonfrag),
                            "uniquefragsa": list(uniquefragsa),
                            "uniquefragsq": list(uniquefragsq),
                        }
                    )
            queryfraginfo = queryfraginfoa[expand]

    return fragcompare


def summarizefrags(row):
    fragcompare = copy.deepcopy(row.Fragcompare)

    commonfrags1 = [
        comfrag
        for query in fragcompare["fraginfo1"]
        for frag in fragcompare["fraginfo1"][query]
        for comfrag in fragcompare["fraginfo1"][query][frag]["commonfrags"]
    ]
    uniquefragsa1 = [
        uniquefraga
        for query in fragcompare["fraginfo1"]
        for frag in fragcompare["fraginfo1"][query]
        for uniquefraga in fragcompare["fraginfo1"][query][frag]["uniquefragsa"]
    ]
    uniquefragsq1 = [
        uniquefragq
        for query in fragcompare["fraginfo1"]
        for frag in fragcompare["fraginfo1"][query]
        for uniquefragq in fragcompare["fraginfo1"][query][frag]["uniquefragsq"]
    ]

    commonfrags2 = [
        comfrag
        for query in fragcompare["fraginfo2"]
        for frag in fragcompare["fraginfo2"][query]
        for comfrag in fragcompare["fraginfo2"][query][frag]["commonfrags"]
    ]
    uniquefragsa2 = [
        uniquefraga
        for query in fragcompare["fraginfo2"]
        for frag in fragcompare["fraginfo2"][query]
        for uniquefraga in fragcompare["fraginfo2"][query][frag]["uniquefragsa"]
    ]
    uniquefragsq2 = [
        uniquefragq
        for query in fragcompare["fraginfo2"]
        for frag in fragcompare["fraginfo2"][query]
        for uniquefragq in fragcompare["fraginfo2"][query][frag]["uniquefragsq"]
    ]

    commonfrags3 = [
        comfrag
        for query in fragcompare["fraginfo3"]
        for frag in fragcompare["fraginfo3"][query]
        for comfrag in fragcompare["fraginfo3"][query][frag]["commonfrags"]
    ]
    uniquefragsa3 = [
        uniquefraga
        for query in fragcompare["fraginfo3"]
        for frag in fragcompare["fraginfo3"][query]
        for uniquefraga in fragcompare["fraginfo3"][query][frag]["uniquefragsa"]
    ]
    uniquefragsq3 = [
        uniquefragq
        for query in fragcompare["fraginfo3"]
        for frag in fragcompare["fraginfo3"][query]
        for uniquefragq in fragcompare["fraginfo3"][query][frag]["uniquefragsq"]
    ]

    return (
        commonfrags1,
        uniquefragsa1,
        uniquefragsq1,
        commonfrags2,
        uniquefragsa2,
        uniquefragsq2,
        commonfrags3,
        uniquefragsa3,
        uniquefragsq3,
    )


### Ranking ###


def rank_impuritiest(
    summary,
    mainrxn,
    impfinalfilt,
    trange=None,
    relquantile=0.75,
    uppertquantile=0.9,
    lowertquantile=0.1,
    process="explode",
):
    #     breakpoint()
    summarytglobal = []
    summaryt = {}
    mainframe = copy.deepcopy(impfinalfilt.loc[impfinalfilt.impurityrxn == mainrxn])
    mainframe2 = mainframe.loc[mainframe.Temperature.astype(bool)]
    mainframe3 = mainframe2.explode("Temperature")
    mainframe3["Temperature"] = mainframe3["Temperature"].apply(pd.to_numeric)
    print(
        "Only "
        + str(round((len(mainframe2.index) / len(mainframe.index)) * 100, -1))
        + "% have temperature records (main reaction)"
    )
    maxtemp = int(mainframe3.Temperature.max())
    mintemp = int(mainframe3.Temperature.min())
    binwidth = 10
    refbins = list(range(mintemp, maxtemp + binwidth, binwidth))
    fig = constructhist(mainframe3.Temperature, refbins, "Main product range")
    barcount, barrel = constructbar(mainframe3, refbins, "Main product range")
    mainframe = mainframe.sort_values(by="Relevance_morgan", ascending=False)
    mainframe2 = mainframe2.sort_values(by="Relevance_morgan", ascending=False)
    mainframe3 = mainframe3.sort_values(by="Relevance_morgan", ascending=False)
    mainframe["Relevance_morgan"] = mainframe["Relevance_morgan"].apply(
        lambda x: round(x, 2)
    )
    mainframe2["Relevance_morgan"] = mainframe2["Relevance_morgan"].apply(
        lambda x: round(x, 2)
    )
    mainframe3["Relevance_morgan"] = mainframe3["Relevance_morgan"].apply(
        lambda x: round(x, 2)
    )
    if process == "explode":
        thresh = mainframe3.Relevance_morgan.quantile(relquantile)
    else:
        thresh = mainframe2.Relevance_morgan.quantile(relquantile)
    lowerbound = mainframe3.loc[
        mainframe3.Relevance_morgan >= thresh
    ].Temperature.quantile(lowertquantile)
    lowerbound = float(floor(lowerbound))
    upperbound = mainframe3.loc[
        mainframe3.Relevance_morgan >= thresh
    ].Temperature.quantile(uppertquantile)
    upperbound = float(ceil(upperbound))
    t_range = [lowerbound, upperbound]
    #     t_range=[round(lowerbound,2),round(upperbound,2)]
    if trange is not None:
        lowerbound = min(trange)
        upperbound = max(trange)
    #     mainframe4=mainframe3.loc[(mainframe3.Temperature<=upperbound) & (mainframe3.Temperature>=lowerbound) & (mainframe3.Relevance_morgan>=thresh)].sort_values(by='Relevance_morgan',ascending=False)
    mainframe4 = mainframe3.loc[
        (mainframe3.Temperature <= upperbound) & (mainframe3.Temperature >= lowerbound)
    ].sort_values(by="Relevance_morgan", ascending=False)
    max_relevance = mainframe4.Relevance_morgan.max()
    rejrxns = mainframe.loc[
        (~mainframe.Temperature.astype(bool))
        & (mainframe.Relevance_morgan > max_relevance)
    ]
    rejrxns.sort_values(by="Relevance_morgan", ascending=False)
    summaryt.update(
        {
            "rxn": mainrxn,
            "querycompds": mainframe3.iloc[0].querycompds,
            "products": mainframe3.iloc[0].impurities,
            "Frame": mainframe,
            "Hits_old": len(mainframe.index),
            "Max relevance_old": mainframe.Relevance_morgan.max(),
            "Frame2": mainframe2,
            "tindicated (%)": round(
                (len(mainframe2.index) / len(mainframe.index)) * 100, -1
            ),
            "Hits_tonly": len(mainframe2.index),
            "Max relevance_tonly": mainframe2.Relevance_morgan.max(),
            "Frame3": mainframe3,
            "thistogram": fig,
            "tbarfreq": barcount,
            "tbarrel": barrel,
            "Hits_tonlyexp": len(mainframe3.index),
            "Relevance threshold": thresh,
            "t_range": t_range,
            "Frame4": mainframe4,
            "Max relevance_tfiltered": max_relevance,
            "Hits_tfiltered": len(mainframe4.index),
            "Reactions_missing": rejrxns,
            "Max relevance_missing": rejrxns.Relevance_morgan.max(),
            "Hits_missing": len(rejrxns.index),
        }
    )
    summarytglobal += [summaryt]
    for imprxn in summary.index.get_level_values(0):
        #         print(imprxn)
        #         breakpoint()
        summaryt = {}
        if imprxn == mainrxn:
            continue
        impframe = copy.deepcopy(impfinalfilt.loc[impfinalfilt.impurityrxn == imprxn])
        impframe2 = impframe.loc[impframe.Temperature.astype(bool)]
        impframe3 = impframe2.explode("Temperature")
        impframe3["Temperature"] = impframe3["Temperature"].apply(pd.to_numeric)
        tindicated = round((len(impframe2.index) / len(impframe.index)) * 100, -1)
        print("Only " + str(tindicated) + "% have temperature records (impurity)")
        if tindicated != 0:  # No temperature entries
            maxtempi = int(impframe3.Temperature.max())
            mintempi = int(impframe3.Temperature.min())
            upperbins = []
            lowerbins = []
            if mintempi < mintemp:
                lowerbins = list(range(mintempi, mintemp, binwidth))
            if maxtempi > maxtemp:
                upperbins = list(
                    range(maxtemp + binwidth, maxtempi + binwidth, binwidth)
                )
            bins = lowerbins + refbins + upperbins
            fig = constructhist(
                mainframe3.Temperature,
                bins,
                "Main product range",
                col2=impframe3.Temperature,
                label2="Impurity range",
            )
            barcount, barrel = constructbar(
                mainframe3,
                bins,
                "Main product range",
                impframe=impframe3,
                label2="Impurity range",
            )
        else:
            fig = None
            barcount = None
            barrel = None
        impframe = impframe.sort_values(by="Relevance_morgan", ascending=False)
        impframe2 = impframe2.sort_values(by="Relevance_morgan", ascending=False)
        impframe3 = impframe3.sort_values(by="Relevance_morgan", ascending=False)
        impframe["Relevance_morgan"] = impframe["Relevance_morgan"].apply(
            lambda x: round(x, 2)
        )
        impframe2["Relevance_morgan"] = impframe2["Relevance_morgan"].apply(
            lambda x: round(x, 2)
        )
        impframe3["Relevance_morgan"] = impframe3["Relevance_morgan"].apply(
            lambda x: round(x, 2)
        )
        if process == "explode":
            thresh = impframe3.Relevance_morgan.quantile(relquantile)
        else:
            thresh = impframe2.Relevance_morgan.quantile(relquantile)
        lowerboundi = impframe3.loc[
            impframe3.Relevance_morgan >= thresh
        ].Temperature.quantile(lowertquantile)
        upperboundi = impframe3.loc[
            impframe3.Relevance_morgan >= thresh
        ].Temperature.quantile(uppertquantile)
        try:
            t_rangei = [float(floor(lowerboundi)), float(ceil(upperboundi))]
        except Exception:
            t_rangei = [lowerboundi, upperboundi]
        #         impframe4=impframe3.loc[(impframe3.Temperature<=upperbound) & (impframe3.Temperature>=lowerbound) & (impframe3.Relevance_morgan>=thresh)].sort_values(by='Relevance_morgan',ascending=False)
        impframe4 = impframe3.loc[
            (impframe3.Temperature <= upperbound)
            & (impframe3.Temperature >= lowerbound)
        ].sort_values(by="Relevance_morgan", ascending=False)
        max_relevance = impframe4.Relevance_morgan.max()
        rejrxns = impframe.loc[
            (~impframe.Temperature.astype(bool))
            & (impframe.Relevance_morgan > max_relevance)
        ]
        rejrxns.sort_values(by="Relevance_morgan", ascending=False)
        summaryt.update(
            {
                "rxn": imprxn,
                "querycompds": impframe.iloc[0].querycompds,
                "products": impframe.iloc[0].impurities,
                "Frame": impframe,
                "Hits_old": len(impframe.index),
                "Max relevance_old": impframe.Relevance_morgan.max(),
                "Frame2": impframe2,
                "tindicated (%)": tindicated,
                "Hits_tonly": len(impframe2.index),
                "Max relevance_tonly": impframe2.Relevance_morgan.max(),
                "Frame3": impframe3,
                "thistogram": fig,
                "tbarfreq": barcount,
                "tbarrel": barrel,
                "Hits_tonlyexp": len(impframe3.index),
                "Relevance threshold": thresh,
                "t_range": t_rangei,
                "Frame4": impframe4,
                "Max relevance_tfiltered": max_relevance,
                "Hits_tfiltered": len(impframe4.index),
                "Reactions_missing": rejrxns,
                "Max relevance_missing": rejrxns.Relevance_morgan.max(),
                "Hits_missing": len(rejrxns.index),
            }
        )
        summarytglobal += [summaryt]
        plt.close("all")
    return summarytglobal


def constructhist(col, bins, label, col2=None, label2=None):
    fig = plt.figure()
    hist_1, bins1, p = plt.hist(
        col,
        bins=bins,
        weights=np.zeros_like(col) + 1.0 / col.size,
        label=label,
        alpha=0.5,
    )
    if col2 is not None:
        hist_2, bins2, p2 = plt.hist(
            col2,
            bins=bins,
            weights=np.zeros_like(col2) + 1.0 / col2.size,
            label=label2,
            alpha=0.5,
        )
    plt.xlabel("Temperature ($^\circ$C)")
    plt.ylabel("Relative Frequency")
    plt.title("Histogram of temperatures")
    plt.legend(loc="best")
    return fig


def constructbar(mainframe, bins, label, impframe=None, label2=None, binwidth=10):
    relm = []
    countm = []
    plotbins = []
    if impframe is not None:
        reli = []
        counti = []
    for refb in bins:
        plotbins += [str(refb) + "\n" + " to " + "\n" + str(refb + binwidth)]
        maxrelm = mainframe.loc[
            (mainframe.Temperature < refb + binwidth) & (mainframe.Temperature >= refb)
        ].Relevance_morgan.max()
        totcountm = mainframe.loc[
            (mainframe.Temperature < refb + binwidth) & (mainframe.Temperature >= refb)
        ].size
        if np.isnan(maxrelm):
            maxrelm = 0
        relm += [maxrelm]
        countm += [totcountm]
        if impframe is not None:
            maxreli = impframe.loc[
                (impframe.Temperature < refb + binwidth)
                & (impframe.Temperature >= refb)
            ].Relevance_morgan.max()
            totcounti = impframe.loc[
                (impframe.Temperature < refb + binwidth)
                & (impframe.Temperature >= refb)
            ].size
            if np.isnan(maxreli):
                maxreli = 0
            reli += [maxreli]
            counti += [totcounti]

    fig_count = plt.figure(figsize=(25, 10))
    plt.rc("font", size=14)
    plt.bar(plotbins, countm, label="Main product range", alpha=0.5)
    if impframe is not None:
        plt.bar(plotbins, counti, label="Impurity range", alpha=0.5)
    plt.xlabel("Temperature range ($^\circ$C)")
    plt.ylabel("Frequency")
    plt.title("Temperature bar chart (frequency)")
    plt.legend(loc="best")

    fig_rel = plt.figure(figsize=(25, 10))
    plt.rc("font", size=14)
    plt.bar(plotbins, relm, label="Main product range", alpha=0.5)
    if impframe is not None:
        plt.bar(plotbins, reli, label="Impurity range", alpha=0.5)
    plt.xlabel("Temperature range ($^\circ$C)")
    plt.ylabel("Max. relevance")
    plt.title("Temperature bar chart (relevance)")
    plt.legend(loc="best")
    return fig_count, fig_rel
