from collections import Counter, OrderedDict
import copy
import itertools
from math import ceil
from typing import Dict, List

import modin.pandas as mpd
import numpy as np
import pandas as pd
from rdkit import Chem  # Importing RDKit
from rdkit.Chem import (
    rdChemReactions,
    rdPartialCharges,
    rdmolops,
)  # Reaction processing

from AnalgCompds import getCarrierFrags0
from AnalgRxns import getspecdat_rxn
from BalanceRxns import balancerxn, checkrxn, maprxn
from GenTempl import clear_isotope, gen_template_row
from helpCompound import hc_Dict
from MainFunctions import (
    drawReaction,
    getcompdict,
    getfragments,
    initray,
    mol_with_atom_index,
    molfromsmiles,
)
from MapRxns import checkrxnrow, update_matches  # , updaterxns_
from RxnCenter import getrxncenterrow, validrxncenterrow
from IPython.display import display

import ray

IP = {
    "reagents": [],  # Only if one reaction is inputted
    "solvents": [],  # Only if one reaction is inputted
    "coefflim": 6,  # Maximum tolerable stoichiometric coefficient
    "usemapper": True,  # If mapper needs to be used
    "addrctonly": False,  # If only reactants should be included for balancing
    "ignoreH": False,  # If hydrogens are to be ignored when balancing
    "hc_prod": hc_Dict,  # Help compound dictionary,
    "coefflim": 6,
    "hc_react": None,
    "first": True,
    "ncpus": 1,
    "restart": True,
    "shutdown_after": False,
}  # Help reactant dictionary


def balance_rxn_uspto_df(df: pd.DataFrame, IP=IP, **kwargs):
    if kwargs:
        IP = {**IP, **kwargs}
    if IP["ncpus"] > 1:
        initray(num_cpus=IP["ncpus"], restart=IP["restart"])
        dfdis = mpd.DataFrame(df)
    else:
        dfdis = df
    dfbal = dfdis.apply(
        balance_rxn_uspto_row,
        IP=IP,
        axis=1,
        result_type="reduce",
    )
    dfser = pd.Series(
        data=dfbal.values, index=dfbal.index
    )  # Optional convert modin back to pandas
    dfbal0 = pd.DataFrame(
        data=dfser.tolist(),
        index=dfser.index,
        columns=[
            "rxnsmiles0",
            "Rdata",
            "Pdata",
            "balrxnsmiles",
            "msg",
            "LHSids",
            "RHSids",
            "hcrct",
            "hcprod",
            "LHSdata",
            "RHSdata",
            "Rgtdata",
            "Solvdata",
            "mappedrxn",
            "conf",
            "msg1",
        ],
    )
    if IP["shutdown_after"] and IP["ncpus"] > 1:
        ray.shutdown()
    return dfbal0


def balance_rxn_uspto_row(row: pd.Series, IP=IP, **kwargs):
    if kwargs:
        IP = {**IP, **kwargs}
    (
        rxnsmiles0,
        Rdata,
        Pdata,
        balrxnsmiles,
        msg,
        LHSids,
        RHSids,
        hcrct,
        hcprod,
        LHSdata,
        RHSdata,
        Rgtdata,
        Solvdata,
        mappedrxn,
        conf,
        msg1,
    ) = balance_rxn_dp(rxn=row["uspto_reaction_smiles"], IP=IP)
    return (
        rxnsmiles0,
        Rdata,
        Pdata,
        balrxnsmiles,
        msg,
        LHSids,
        RHSids,
        hcrct,
        hcprod,
        LHSdata,
        RHSdata,
        Rgtdata,
        Solvdata,
        mappedrxn,
        conf,
        msg1,
    )


def parse_uspto(
    rxnsmiles: str,
):  # NOTE: MIXTURES NOT WELL REPRESENTED HERE AS REACTION SMILES IS INSUFFICIENT
    reagents = []
    if ">>" not in rxnsmiles:  # Reagents are present
        rxnsmilesgroup = rxnsmiles.split(">")
        reagents = rxnsmilesgroup[1].split(".")
        rxnsmiles = ">>".join([rxnsmilesgroup[0], rxnsmilesgroup[2]])
    return rxnsmiles, reagents


def balance_rxn_dp(rxn: str = "", IP=IP, rctstorgts=True, **kwargs):
    """_summary_

    Args:
        rxn (str): _description_

    Returns:
        _type_: _description_
    """
    if kwargs:
        IP = {**IP, **kwargs}
    if IP["hc_prod"] is None:
        IP["hc_prod"] = {}
    Rdata = {}
    Rgtdata = {}
    Solvdata = {}
    Pdata = {}
    dicts = [
        "Rdata",
        "Pdata",
        "Rgtdata",
        "Solvdata",
    ]  # Allowable inputs, either {} or {..} output of getcompdict
    try:  # First try with user inputs
        if "Rgtdata" in IP:
            Rgtdata = IP["Rgtdata"]
            dicts.remove("Rgtdata")
        if "Solvdata" in IP:
            Solvdata = IP["Solvdata"]
            dicts.remove("Solvdata")
        if "Rdata" in IP:
            Rdata = IP["Rdata"]
            dicts.remove("Rdata")
        if "Pdata" in IP:
            Pdata = IP["Pdata"]
            dicts.remove("Pdata")
        if dicts:
            if rxn:  # Reaction string is specified
                if ">>" not in rxn:  # with reagents eg. USPTO)
                    rxn, reagents = parse_uspto(rxn)
                    IP["reagents"] = reagents

                output = getspecdat_rxn(
                    rxn, reagents=IP["reagents"], solvents=IP["solvents"], dicts=dicts
                )
                if "Rdata" in dicts:
                    Rdata = output[0]
                if "Pdata" in dicts:
                    Pdata = output[1]
                if "Rgtdata" in dicts:
                    Rgtdata = output[2]
                if "Solvdata" in dicts:
                    Solvdata = output[3]

    except Exception as e:
        msg = "Invalid. Species missing. " + str(e)
        return (
            "Error",
            "Error",
            "Error",
            "Error",
            msg,
            [],
            [],
            [],
            [],
            "Error",
            "Error",
            {},
            {},
            "Error",
            "Error",
            msg,
        )
    # IP["addedspecies"] = [i for i in Rdata]
    rxnsmiles0 = ">>".join(
        [
            getfragments(
                [Rdata[r]["smiles"] for r in Rdata for _ in range(Rdata[r]["count"])],
                smiles=True,
            ),
            getfragments(
                [Pdata[p]["smiles"] for p in Pdata for _ in range(Pdata[p]["count"])],
                smiles=True,
            ),
        ]
    )
    # New
    if "msg" in IP:
        if IP["msg"] and "With hydrogen carriers" in IP["msg"]:
            IP["msg"] = "With hydrogen carriers: " + ",".join(
                [
                    hcarrier
                    for hcarrier in IP["msg"]
                    .split("With hydrogen carriers: ")[1]
                    .split(", ")[0]
                    .split(",")
                    if hcarrier.isdigit()
                ]
            )
        else:
            IP["msg"] = ""
    IP["rxnsmiles0"] = rxnsmiles0
    input = {
        key: IP[key]
        for key in [
            "rxnsmiles0",
            "first",
            "usemapper",
            "addedspecies",
            "hc_prod",
            "hc_react",
            "coefflim",
            "addrctonly",
            "ignoreH",
            "mandrcts",
            "msg",
        ]
        if key in IP
    }

    (
        rxnsmiles0,
        balrxnsmiles,
        msg,
        LHSids,
        RHSids,
        hcrct,
        hcprod,
        LHSdata,
        RHSdata,
    ) = balancerxn(
        Rdata,
        Pdata,
        Rgtdata=Rgtdata,
        Solvdata=Solvdata,
        **input,
    )
    # print(msg)
    # print(balrxnsmiles)
    # print(Rdata)
    # print(LHSdata)
    mappedrxn_ = maprxn([balrxnsmiles])[0]
    if mappedrxn_ != "Error":
        mappedrxn = mappedrxn_.get("mapped_rxn")
        conf = mappedrxn_.get("confidence")
    else:
        mappedrxn = "Error"
        conf = 0
    if mappedrxn != "Error":
        rseries = pd.DataFrame(
            [
                {
                    "mapped_rxn": mappedrxn,
                    "LHSdata": LHSdata,
                    "RHSdata": RHSdata,
                    "msg": msg,
                }
            ]
        ).iloc[0]
        LHSdata, RHSdata, msg1 = checkrxnrow(rseries)
        if ("Unmapped" in msg1 or "unmapped" in msg1) and (
            "Smiles discrepancy" not in msg1
        ):  # Need to skip if mandatory reactants carries over
            if (
                rctstorgts and "unmapped" in msg1
            ):  # Adding unmapped mandatory reactants as reagents
                Rgtdata.update(
                    {
                        int(ID): copy.deepcopy(Rdata[int(ID)])
                        for ID in msg1.split("unmapped from LHS: ")[1]
                        .split(", ")[0]
                        .split(",")
                        if ID.isdigit()
                    }
                )

            inputdict = {
                "Rdata": Rdata,
                "LHSdata": LHSdata,
                "RHSdata": RHSdata,
                "msg": msg,
                "msg1": msg1,
                "hcprod": hcprod,
                "hcrct": hcrct,
                "Rgtdata": Rgtdata,
                "Solvdata": Solvdata,
            }
            if (
                "mandrcts" in input and input["mandrcts"]
            ):  # External mandatory reactants supplied
                inputdict["Rdata"] = input["mandrcts"]
            # print(inputdict)

            (
                mappedrxn,
                conf,
                balrxnsmiles,
                msg,
                LHSids,
                RHSids,
                hcrct,
                hcprod,
                LHSdata,
                RHSdata,
                msg1,
            ) = updaterxns_(pd.DataFrame([inputdict]).iloc[0], hc_prod=IP["hc_prod"])
    else:
        msg1 = "Mapping error"
    LHSids = [ID for ID in LHSdata for _ in range(LHSdata[ID]["count"])]
    RHSids = [ID for ID in RHSdata for _ in range(RHSdata[ID]["count"])]
    return (
        rxnsmiles0,
        Rdata,
        Pdata,
        balrxnsmiles,
        msg,
        LHSids,
        RHSids,
        hcrct,
        hcprod,
        LHSdata,
        RHSdata,
        Rgtdata,
        Solvdata,
        mappedrxn,
        conf,
        msg1,
    )


# balance_rxn_dp('CCCCCCc1ccc(C(=O)Cl)cc1.CCCCc1ccc(C#Cc2ccc(CNc3ccc4c(c3)C(=O)OC(C)(C)O4)cc2)cc1.Cl>>CCCCCCc1ccc(C(=O)N(Cc2ccc(C#Cc3ccc(CCCC)cc3)cc2)c2ccc3c(c2)C(=O)OC(C)(C)O3)cc1')


def rxn_center_dp(mappedrxn: str, LHSdata: Dict, RHSdata: Dict, expand: int = 1):
    # print(LHSdata)
    # print(RHSdata)
    for rctid in LHSdata:
        userinput = LHSdata[rctid]["smiles"]
        fraglist = getCarrierFrags0(userinput, resFormat="smiles", expand=expand)
        fragloc = {}
        nofg = set()
        for idx, cleanmol in enumerate(LHSdata[rctid]["cleanmol"]):
            # cleanmol = Chem.AddHs(cleanmol)
            for frag in fraglist:
                fragloc, nofg = update_matches(
                    Chem.AddHs(Chem.RemoveAllHs(cleanmol)),
                    frag,
                    fragloc=fragloc,
                    nofg=nofg,
                    idx=idx,
                    rctid=rctid,
                )
        LHSdata[rctid]["fragloc"] = fragloc
    specmap, rnbmap, rxncentermapnum, msg = getrxncenterrow(
        pd.DataFrame(
            [{"mapped_rxn": mappedrxn, "LHSdata": LHSdata, "RHSdata": RHSdata}]
        ).iloc[0],
    )
    if msg:
        LHSdata, msg, outfrag, outfg, outneighbor, unusedanalogue = validrxncenterrow(
            pd.DataFrame(
                [
                    {
                        "specmap": specmap,
                        "rxncentermapnum": rxncentermapnum,
                        "LHSdata": LHSdata,
                        "rnbmap": rnbmap,
                    }
                ]
            ).iloc[0]
        )
    else:
        outfrag = "Error"
        outfg = "Error"
        outneighbor = "Error"
        unusedanalogue = "Error"
    return (
        specmap,
        rnbmap,
        rxncentermapnum,
        LHSdata,
        msg,
        outfrag,
        outfg,
        outneighbor,
        unusedanalogue,
    )


def gen_template_ip(
    LHSdata: Dict,
    RHSdata: Dict,
    specmap: Dict,
    outfrag: Dict = {},
    rnbmap: Dict = {},
    unusedanalogue: List = [],
    specificity="loose",
    processall=True,
):
    template, LHSdata, RHSdata, msg4, farfg, unusedprod = gen_template_row(
        pd.DataFrame(
            [
                {
                    "LHSdata": LHSdata,
                    "RHSdata": RHSdata,
                    "specmap": specmap,
                    "outfrag": outfrag,
                    "rnbmap": rnbmap,
                    "unusedanalogue": unusedanalogue,
                }
            ]
        ).iloc[0],
        specificity=specificity,
        processall=processall,
    )
    return template, LHSdata, RHSdata, msg4, farfg, unusedprod


def updaterxns_(row, hc_prod={}):
    """
    Updates reactions if there are unmapped species and balances if there are changes (optional). Assumes both
    balancerxn, maprxn and checkrxn have all been called already

    rctstorgts: Move unmapped mandatory reactants to reagents

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
    storemsg = ""
    i = 0
    if "Smiles discrepancy" in msg1:
        raise Exception("Smiles discrepancy")
    while (
        "Unmapped" in msg1 or "unmapped" in msg1
    ) or i == 0:  # Unmapped species exist not reflected
        if "from RHS" or "from LHS" in msg1:  # Mandatory products/reactants unmapped
            if "from LHS" in msg1 and i != 0:  # To avoid infinite loop
                break
            storemsg = msg1
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
            elif msg1 != storemsg:
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


# def comparetemplate(reftemplate: str, template: str):
#     reftemplrxn = rdChemReactions.ReactionFromSmarts(template)
#     templrxn = rdChemReactions.ReactionFromSmarts(template)


def visrxn(row):
    print(row.rxnsmiles0)
    display(drawReaction(row.rxnsmiles0))
    print(row.mappedrxn)
    display(drawReaction(row.mappedrxn))


def parserow(row):
    rxnsmiles0 = row.rxnsmiles0
    Rdata = row.Rdata
    Pdata = row.Pdata
    balrxnsmiles = row.balrxnsmiles
    msg = row.msg
    LHSids = row.LHSids
    RHSids = row.RHSids
    hcrct = row.hcrct
    hcprod = row.hcprod
    LHSdata = row.LHSdata
    RHSdata = row.RHSdata
    Rgtdata = row.Rgtdata
    Solvdata = row.Solvdata
    mappedrxn = row.mappedrxn
    conf = row.conf
    msg1 = row.msg1
    return (
        rxnsmiles0,
        Rdata,
        Pdata,
        balrxnsmiles,
        msg,
        LHSids,
        RHSids,
        hcrct,
        hcprod,
        LHSdata,
        RHSdata,
        Rgtdata,
        Solvdata,
        mappedrxn,
        conf,
        msg1,
    )


def fullsmiles(rxnsmiles, Rgtdata={}, Solvdata={}):
    Rgtdata = {**Rgtdata, **Solvdata}
    if Rgtdata:
        rxnlist = rxnsmiles.split(">>")
        fullrxn = (
            rxnlist[0]
            + ">"
            + ".".join(
                [
                    Rgtdata[ID]["smiles"]
                    for ID in Rgtdata
                    for _ in range(Rgtdata[ID]["count"])
                    if Rgtdata[ID]["smiles"] not in rxnlist[0].split(".")
                ]
            )
            + ">"
            + rxnlist[1]
        )
    return fullrxn


def parsecandidate(candidate, neutralizecharge=False):
    """Parses candidate, adds hydrogen if possible and sanitizes product

    Args:
        candidate (str): Candidate SMILES

    Returns:
        str: Sanitized SMILES
    """
    candidatemol = Chem.MolFromSmarts(candidate)
    candidatemol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(candidatemol)
    if neutralizecharge and rdmolops.GetFormalCharge(candidatemol) != 0:
        candidatemol = neutralize_atoms(candidatemol)
    candidate = Chem.MolToSmiles(candidatemol)
    return candidate


def resolvechargebal(
    LHSdata, RHSdata, refdict, Rdata={}, **kwargs
):  # Only resolves LHS, refdict is Rdata for USPTO and Rgtdata+Rdata for Reaxys
    LHSdata2 = copy.deepcopy(LHSdata)
    RHSdata2 = copy.deepcopy(RHSdata)
    Rcharge = sum(
        [
            LHSdata2[ID]["charge"]
            for ID in LHSdata2
            for _ in range(LHSdata2[ID]["count"])
        ]
    )
    Pcharge = sum(
        [
            RHSdata2[ID]["charge"]
            for ID in RHSdata2
            for _ in range(RHSdata2[ID]["count"])
        ]
    )
    if Rcharge != Pcharge:
        chargedict = {}
        for ccandi in refdict:
            charge = refdict[ccandi]["charge"] * refdict[ccandi]["count"]
            if (Rcharge < 0 and charge > 0) or (Rcharge > 0 and charge < 0):
                if charge not in chargedict:
                    chargedict[charge] = [ccandi]
                else:
                    chargedict[charge].append(ccandi)
        chosen_comb = []
        combfound = False
        if len(chargedict) == 1:
            charge = list(chargedict.keys())[0]
            mult = abs(Rcharge / charge)
            if mult.is_integer:
                combination = [charge for _ in range(int(mult))]
                counter_ = Counter(combination)
                chosen_comb = counter_
                combfound = True
        else:
            filteredclist = []
            if Rdata:
                for ID in Rdata:
                    if Rdata[ID]["charge"] != 0 and ID in LHSdata2:
                        oldcount = Rdata[ID]["count"]
                        newcount = LHSdata2[ID]["count"]
                        if newcount > oldcount:
                            filteredclist += [-1 * Rdata[ID]["charge"]]
            if filteredclist:
                chargedict0 = copy.deepcopy(chargedict)
                chargedict = {k: v for k, v in chargedict.items() if k in filteredclist}
                if not chargedict:
                    chargedict = chargedict0
            minidx = ceil(abs(Rcharge) / max(chargedict.keys()))
            maxidx = ceil(abs(Rcharge) / min(chargedict.keys()))
            combinations = []
            combinations_counter = []
            for i in range(minidx, maxidx + 1):
                combinations = [
                    c
                    for c in itertools.combinations_with_replacement(
                        list(chargedict.keys()), i
                    )
                    if sum(c) == (Rcharge * -1)
                ]
                #             print(combinations)
                if combinations:
                    for combination in combinations:
                        counter_ = Counter(combination)
                        combinations_counter += [counter_]
                        if all([counter_[key] == 1 for key in counter_]):
                            print("Combination found")
                            combfound = True
                            chosen_comb = counter_
                            break
                if combfound:
                    break
        if not combfound:
            prodarray = np.array(
                [
                    np.prod(np.array(list(counter_.values())))
                    for counter_ in combinations_counter
                ]
            )
            idxmin = np.where(prodarray == prodarray.min())
            combinations_counter = [
                combinations_counter[idxmin_] for idxmin_ in idxmin[0]
            ]
            if len(combinations_counter) > 1:
                lenarray = np.array(
                    [len(counter_) for counter_ in combinations_counter]
                )
                idxmin2 = np.where(lenarray == lenarray.min())
                combinations_counter = [
                    combinations_counter[idxmin_] for idxmin_ in idxmin2[0]
                ]
            if len(combinations_counter) == 1:
                print("Combination found")
                combfound = True
                chosen_comb = combinations_counter[idxmin]
        if chosen_comb:
            for charge in chosen_comb:
                ccandi = chargedict[charge]
                if len(ccandi) > 1:
                    print("Too many charge candidates")
                    raise Exception("Too many charge candidates")
                else:
                    ccandi = ccandi[0]
                    LHSdata2.update({ccandi: refdict[ccandi]})
                    LHSdata2[ccandi]["count"] = (
                        chosen_comb[charge] * LHSdata2[ccandi]["count"]
                    )
            LHSsmiles0 = getfragments_dict(LHSdata2, usemapped=False)
            RHSsmiles0 = getfragments_dict(RHSdata2, usemapped=False)
            LHSsmiles1 = getfragments_dict(LHSdata2)
            RHSsmiles1 = getfragments_dict(RHSdata2)
            rxnsmiles1 = LHSsmiles1 + ">>" + RHSsmiles1
            rxn = rdChemReactions.ReactionFromSmarts(rxnsmiles1)
            rxnsmiles0 = LHSsmiles0 + ">>" + RHSsmiles0
            return cleavagebal(
                LHSdata2,
                RHSdata2,
                rxnsmiles0,
                rxn,
                resolvechargebal=False,
                Rdata=Rdata,
                hc_prod={},
                neutralizecharge=True,
                **kwargs,
            )
        else:
            raise Exception("No possible candidates available")


def cleavagebal(
    LHSdata,
    RHSdata,
    rxnsmiles,
    rxn,
    msg="",
    Rdata={},
    Rgtdata={},
    Solvdata={},
    resolvecharge=False,
    neutralizecharge=False,
    **kwargs
):
    """_summary_

    Args:
        LHSdata (_type_): _description_
        RHSdata (_type_): _description_
        rxnsmiles (_type_): Unmapped smiles
        rxn (_type_): Mapped smiles
        msg (str, optional): _description_. Defaults to ''.
        Rdata (dict, optional): _description_. Defaults to {}.
        Rgtdata (dict, optional): _description_. Defaults to {}.
        Solvdata (dict, optional): _description_. Defaults to {}.
        resolvechargebal (bool, optional): _description_. Defaults to False.
    """

    finallist = []
    oldfinallist = []
    LHSdata2 = copy.deepcopy(LHSdata)
    RHSdata2 = copy.deepcopy(RHSdata)
    if "LHSids" in kwargs and kwargs["LHSids"]:
        LHSids = kwargs["LHSids"]
    else:
        LHSids = [ID for ID in LHSdata2 for _ in range(LHSdata2[ID]["count"])]
    for j, rct in enumerate(rxn.GetReactants()):
        unmappedatoms = [
            atom.GetIdx()
            for atom in rct.GetAtoms()
            if not atom.HasProp("molAtomMapNumber")
        ]
        if len(unmappedatoms) == len(rct.GetAtoms()):  # Entire reactant is unmapped
            if (
                msg and "With hydrogen carriers:" in msg
            ):  # CRUCIAL, otherwise infinite loop possible
                hcarriers = [
                    int(ID)
                    for ID in msg.split("With hydrogen carriers: ")[1]
                    .split(", ")[0]
                    .split(",")
                    if ID.isdigit()
                ]
                if (
                    LHSids[j] in hcarriers
                ):  # Hydrogen carrier for hydrogen balance..exempt from this function
                    continue
        # print(rxnsmiles)
        z = [unmappedatoms]
        print(z)
        if not any(z):
            print("No candidates present")
            continue
        oldcandidatelist = []
        candidatelist = []
        try:
            for subfrag in z:
                candidate = Chem.MolFragmentToSmiles(rct, atomsToUse=subfrag)
                if "." in candidate:
                    listcandidate = candidate.split(".")
                    oldcandidatelist += listcandidate
                    for candidate_ in listcandidate:
                        candidate = parsecandidate(
                            candidate_, neutralizecharge=neutralizecharge
                        )
                        candidatelist += [candidate]
                else:
                    oldcandidatelist += [candidate]
                    candidate = parsecandidate(
                        candidate, neutralizecharge=neutralizecharge
                    )
                    candidatelist += [candidate]
        except Exception as e:
            print(e)
        finallist += [candidatelist]
        oldfinallist += [oldcandidatelist]
    print(finallist)
    print(oldfinallist)

    # For Reaxys, when species IDs are known
    maxidx = max(RHSdata2.keys())
    masterlist = [
        candidate for candidatelist in finallist for candidate in candidatelist
    ]
    oldmasterlist = [
        candidateold
        for candidatelistold in oldfinallist
        for candidateold in candidatelistold
    ]
    masterlistcounter = Counter(masterlist)
    addedcandidates = []
    idx = maxidx
    for i, candidate in enumerate(masterlist):
        if candidate not in addedcandidates:
            idx += 1
            RHSdata2.update(getcompdict(ID=idx, smiles=candidate))
            RHSdata2[idx]["count"] = masterlistcounter[candidate]
            RHSdata2[idx]["unmappedfrag"] = oldmasterlist[i]
            addedcandidates += [candidate]

    # rxnlist=rxnsmiles.split('>>')
    # rxnlist[1]=rxnlist[1]+'.'+'.'.join([candidate_ for candidate in finallist for candidate_ in candidate])
    # rxnsmiles02='>>'.join(rxnlist)

    (
        rxnsmiles02,
        _,
        _,
        balrxnsmiles2,
        msg0,
        LHSids2,
        RHSids2,
        hcrct2,
        hcprod2,
        LHSdata2,
        RHSdata2,
        _,
        _,
        mappedrxn2,
        conf2,
        msg1_0,
    ) = balance_rxn_dp(
        Rdata=LHSdata2,
        Pdata=RHSdata2,
        Rgtdata=Rgtdata,
        Solvdata=Solvdata,
        ignoreH=True,
        msg=msg,
        hc_prod={},
        first=False,
    )
    print(balrxnsmiles2)
    print(LHSdata2)
    print(msg0)
    print(msg1_0)

    if resolvecharge and "Charge imbalance" in msg0:
        print("Resolving charge imbalance")
        try:
            (
                rxnsmiles02,
                balrxnsmiles2,
                msg0,
                LHSids2,
                RHSids2,
                hcrct2,
                hcprod2,
                LHSdata2,
                RHSdata2,
                _,
                _,
                mappedrxn2,
                conf2,
                msg1_0,
            ) = resolvechargebal(
                LHSdata2, RHSdata2, Rdata, Rdata=Rdata, ignoreH=True, msg=msg0
            )
        except Exception as e:
            # raise Exception(e)
            print(str(e))
            pass

    # Resolving Rgtdata/Solvdata
    if Rgtdata:
        Rgtdata2 = copy.deepcopy(Rgtdata)
        Rgtdata2 = {ID: Rgtdata[ID] for ID in Rgtdata if ID not in LHSids2}
    else:
        Rgtdata2 = {}
    if Solvdata:
        Solvdata2 = copy.deepcopy(Solvdata)
        Solvdata2 = {ID: Solvdata[ID] for ID in Solvdata if ID not in LHSids2}
    else:
        Solvdata2 = {}
    removedmandrcts = []
    # Resolving message
    if Rdata:  # For USPTO if original dictionary is given
        addedspecies = [ID for ID in LHSdata2 if ID not in Rdata]
        removedmandrcts = [ID for ID in Rdata if ID not in LHSdata2]
    else:
        addedspecies = [ID for ID in LHSdata2 if ID in Rgtdata]
    if addedspecies:
        addedstr = ",".join([str(species) for species in set(addedspecies)])
        if "Already balanced" in msg0:
            msg0 = msg0.replace(
                "Already balanced", "Balanced with species: " + addedstr
            )
        elif "Balanced" in msg0:
            msg0 = msg0.replace("Balanced", "Balanced with species: " + addedstr)
        elif " with species" in msg0:
            msg0 = msg0.replace(
                msg0.split(" with species: ")[1].split(", ")[0], addedstr
            )
    # Resolving added RHS species
    hclist = [
        RHSdata2[ID]["formula"] for ID in RHSdata2 if "unmappedfrag" in RHSdata2[ID]
    ]
    if hclist:
        addedstr = ",".join(hclist)
        msg0 = msg0 + ", with cleavage help products: " + addedstr
    # Resolving unmapped mandatory reactants
    if removedmandrcts:
        addedstr = ",".join([str(species) for species in set(removedmandrcts)])
        if "unmapped from LHS" in msg1_0:
            msg1_0 = msg1_0.replace(
                msg1_0.split("unmapped from LHS: ")[1].split(", ")[0], addedstr
            )
        else:
            msg1_0 = "Mandatory species unmapped from LHS: " + addedstr + ", " + msg1_0
    return (
        rxnsmiles02,
        balrxnsmiles2,
        msg0,
        LHSids2,
        RHSids2,
        hcrct2,
        hcprod2,
        LHSdata2,
        RHSdata2,
        Rgtdata2,
        Solvdata2,
        mappedrxn2,
        conf2,
        msg1_0,
    )

    # if 'hydrogen' or 'Hydrogen' in msg0: #Hydrogen surplus on LHS (not sure which species would donate)..try combination


#     if len(candidatelist)==2:

# def combinebal
def getfragments_dict(refdict, usemapped=True):
    if usemapped:
        smiles = []
        for s in refdict:
            if "mappedsmiles" in refdict[s]:
                smiles += refdict[s]["mappedsmiles"]
            else:
                smiles += [refdict[s]["smiles"] for _ in range(refdict[s]["count"])]
        return ".".join(smiles)
    else:
        return getfragments(
            [refdict[s]["smiles"] for s in refdict for _ in range(refdict[s]["count"])],
            smiles=True,
        )


def cleavagebal0(
    LHSdata,
    RHSdata,
    balrxnsmiles,
    mappedsmiles,
    msg="",
    Rdata={},
    Rgtdata={},
    Solvdata={},
    resolvecharge=False,
    **kwargs
):
    """_summary_

    Args:
        LHSdata (_type_): _description_
        RHSdata (_type_): _description_
        rxnsmiles (_type_): Unmapped smiles
        mappedsmiles (_type_): Mapped smiles
        msg (str, optional): _description_. Defaults to ''.
        Rdata (dict, optional): _description_. Defaults to {}.
        Rgtdata (dict, optional): _description_. Defaults to {}.
        Solvdata (dict, optional): _description_. Defaults to {}.
        resolvechargebal (bool, optional): _description_. Defaults to False.
    """

    finallist = []
    oldfinallist = []
    LHSdata2 = copy.deepcopy(LHSdata)
    RHSdata2 = copy.deepcopy(RHSdata)
    if isinstance(mappedsmiles, str):
        rxn = rdChemReactions.ReactionFromSmarts(mappedsmiles)
    else:
        rxn = mappedsmiles
    if msg:  # CRUCIAL, otherwise infinite loop possible
        hcarriers = [
            int(ID)
            for ID in msg.split("With hydrogen carriers: ")[1].split(", ")[0].split(",")
            if ID.isdigit()
        ]

    for j, rct in enumerate(rxn.GetReactants()):
        unmappedatoms = [
            atom.GetIdx()
            for atom in rct.GetAtoms()
            if not atom.HasProp("molAtomMapNumber")
        ]
        if len(unmappedatoms) == len(rct.GetAtoms()):  # Entire reactant is unmapped
            if msg:  # CRUCIAL, otherwise infinite loop possible
                hcarriers = [
                    int(ID)
                    for ID in msg.split("With hydrogen carriers: ")[1]
                    .split(", ")[0]
                    .split(",")
                    if ID.isdigit()
                ]
                if (
                    LHSids[j] in hcarriers
                ):  # Hydrogen carrier for hydrogen balance..exempt from this function
                    continue
        print(balrxnsmiles)
        z = [unmappedatoms]
        print(z)
        if not any(z):
            print("No candidates present")
            continue
        oldcandidatelist = []
        candidatelist = []
        try:
            for subfrag in z:
                candidate = Chem.MolFragmentToSmiles(rct, atomsToUse=subfrag)
                if "." in candidate:
                    listcandidate = candidate.split(".")
                    oldcandidatelist += listcandidate
                    for candidate_ in listcandidate:
                        candidate = parsecandidate(candidate_)
                        candidatelist += [candidate]
                else:
                    oldcandidatelist += [candidate]
                    candidate = parsecandidate(candidate)
                    candidatelist += [candidate]
        except Exception as e:
            print(e)
        finallist += [candidatelist]
        oldfinallist += [oldcandidatelist]
    print(finallist)
    print(oldfinallist)

    # For Reaxys, when species IDs are known
    maxidx = max(RHSdata2.keys())
    masterlist = [
        candidate for candidatelist in finallist for candidate in candidatelist
    ]
    oldmasterlist = [
        candidateold
        for candidatelistold in oldfinallist
        for candidateold in candidatelistold
    ]
    masterlistcounter = Counter(masterlist)
    addedcandidates = []
    idx = maxidx
    for i, candidate in enumerate(masterlist):
        if candidate not in addedcandidates:
            idx += 1
            RHSdata2.update(getcompdict(ID=idx, smiles=candidate))
            RHSdata2[idx]["count"] = masterlistcounter[candidate]
            RHSdata2[idx]["unmappedfrag"] = oldmasterlist[i]
            addedcandidates += [candidate]

    # rxnlist=rxnsmiles.split('>>')
    # rxnlist[1]=rxnlist[1]+'.'+'.'.join([candidate_ for candidate in finallist for candidate_ in candidate])
    # rxnsmiles02='>>'.join(rxnlist)

    (
        rxnsmiles02,
        _,
        _,
        balrxnsmiles2,
        msg0,
        LHSids2,
        RHSids2,
        hcrct2,
        hcprod2,
        LHSdata2,
        RHSdata2,
        _,
        _,
        mappedrxn2,
        conf2,
        msg1_0,
    ) = balance_rxn_dp(
        Rdata=LHSdata2,
        Pdata=RHSdata2,
        Rgtdata=Rgtdata,
        Solvdata=Solvdata,
        ignoreH=True,
        msg=msg,
        hc_prod={},
        first=False,
    )
    print(balrxnsmiles2)
    print(LHSdata2)
    print(msg0)
    print(msg1_0)

    if resolvecharge and "Charge imbalance" in msg0:
        print("Resolving charge imbalance")
        try:
            (
                rxnsmiles02,
                balrxnsmiles2,
                msg0,
                LHSids2,
                RHSids2,
                hcrct2,
                hcprod2,
                LHSdata2,
                RHSdata2,
                _,
                _,
                mappedrxn2,
                conf2,
                msg1_0,
            ) = resolvechargebal(
                LHSdata2, RHSdata2, Rdata, Rdata=Rdata, ignoreH=True, msg=msg0
            )
        except Exception as e:
            # raise Exception(e)
            print(str(e))
            pass

    # Resolving Rgtdata/Solvdata
    if Rgtdata:
        Rgtdata2 = copy.deepcopy(Rgtdata)
        Rgtdata2 = {ID: Rgtdata[ID] for ID in Rgtdata if ID not in LHSids2}
    else:
        Rgtdata2 = {}
    if Solvdata:
        Solvdata2 = copy.deepcopy(Solvdata)
        Solvdata2 = {ID: Solvdata[ID] for ID in Solvdata if ID not in LHSids2}
    else:
        Solvdata2 = {}
    removedmandrcts = []
    # Resolving message
    if Rdata:  # For USPTO if original dictionary is given
        addedspecies = [ID for ID in LHSdata2 if ID not in Rdata]
        removedmandrcts = [ID for ID in Rdata if ID not in LHSdata2]
    else:
        addedspecies = [ID for ID in LHSdata2 if ID in Rgtdata]
    if addedspecies:
        addedstr = ",".join([str(species) for species in set(addedspecies)])
        if "Already balanced" in msg0:
            msg0 = msg0.replace(
                "Already balanced", "Balanced with species: " + addedstr
            )
        elif "Balanced" in msg0:
            msg0 = msg0.replace("Balanced", "Balanced with species: " + addedstr)
        elif " with species" in msg0:
            msg0 = msg0.replace(
                msg0.split(" with species: ")[1].split(", ")[0], addedstr
            )
    # Resolving added RHS species
    hclist = [
        RHSdata2[ID]["formula"] for ID in RHSdata2 if "unmappedfrag" in RHSdata2[ID]
    ]
    if hclist:
        addedstr = ",".join(hclist)
        msg0 = msg0 + ", with cleavage help products: " + addedstr
    # Resolving unmapped mandatory reactants
    if removedmandrcts:
        addedstr = ",".join([str(species) for species in set(removedmandrcts)])
        if "unmapped from LHS" in msg1_0:
            msg1_0 = msg1_0.replace(
                msg1_0.split("unmapped from LHS: ")[1].split(", ")[0], addedstr
            )
        else:
            msg1_0 = "Mandatory species unmapped from LHS: " + addedstr + ", " + msg1_0
    return (
        rxnsmiles02,
        balrxnsmiles2,
        msg0,
        LHSids2,
        RHSids2,
        hcrct2,
        hcprod2,
        LHSdata2,
        RHSdata2,
        Rgtdata2,
        Solvdata2,
        mappedrxn2,
        conf2,
        msg1_0,
    )


def getunmappedfrags(mappedsmiles, refdict, specid=0, inst=0, hcarriers=[]):
    if isinstance(mappedsmiles, str):
        mappedmol = molfromsmiles(mappedsmiles)
    else:
        mappedmol = mappedsmiles
    mappedmol2 = copy.deepcopy(mappedmol)
    mol_with_atom_index(mappedmol2)
    unmappedatoms = [
        atom.GetIdx()
        for atom in mappedmol.GetAtoms()
        if not atom.HasProp("molAtomMapNumber")
    ]
    mappedatoms = [
        atom.GetIdx()
        for atom in mappedmol.GetAtoms()
        if atom.HasProp("molAtomMapNumber")
    ]
    if len(unmappedatoms) == len(mappedmol.GetAtoms()):
        if (
            specid in hcarriers
        ):  # Hydrogen carrier for hydrogen balance..exempt from this function
            print("Hydrogen carrier detected, skipping")
            return refdict
    fracmapped = len(mappedatoms) / (len(unmappedatoms) + len(mappedatoms))
    _, nbmap = storeatommap(mappedmol, specid=specid, inst=inst, inverse=True)
    unmappedatoms2 = []
    if len(nbmap.keys()) > 1:  # Multiple separate fragments
        for nb in nbmap:  # Quite inefficient
            atoms_ = {nb}
            nbs = {nb}
            while nbs:
                nbs = [
                    nb__.GetIdx()
                    for nb_ in nbs
                    for nb__ in mappedmol.GetAtomWithIdx(nb_).GetNeighbors()
                    if nb__.GetIdx() in unmappedatoms
                    if nb__.GetIdx() not in atoms_
                ]
                atoms_ = atoms_.union(set(nbs))
            unmappedatoms2 += [list(atoms_)]
    else:
        unmappedatoms2 += [unmappedatoms]

    if (
        "unmappedatoms" not in refdict[specid]
        or "mappedatoms" not in refdict[specid]
        or "fracmapped" not in refdict[specid]
    ):
        refdict[specid].update(
            {
                "unmappedatoms": {inst: unmappedatoms2},
                "mappedatoms": {inst: mappedatoms},
                "fracmapped": {inst: fracmapped},
            }
        )
    elif inst not in refdict[specid]["unmappedatoms"]:
        refdict[specid]["unmappedatoms"].update({inst: unmappedatoms2})
        refdict[specid]["mappedatoms"].update({inst: mappedatoms})
        refdict[specid]["fracmapped"].update({inst: fracmapped})
    if not any(unmappedatoms2):
        print("No candidates found. All atoms mapped.")
        return refdict
    candilist_simple = []  # Generic representation of the fragment
    candidatelist = []
    rdPartialCharges.ComputeGasteigerCharges(
        mappedmol
    )  # Compute partial charges in the molecule
    try:
        for nb, subfrag in zip(nbmap, unmappedatoms2):
            candidate_simple = Chem.MolFragmentToSmiles(
                mappedmol, atomsToUse=subfrag
            )  # Generates a very generic representation of the fragment
            candidate_strict, atomsymbols, _ = gen_template_fragment(
                set(subfrag),
                mappedmol2,
                specificity="strict",
                seoarateAllHs=True,
                atomsymbols={},
            )  # Generates a strict template fragment, assume mapped mol same as reaction
            charge_m = mappedmol.GetAtomWithIdx(nbmap[nb][2]).GetPropsAsDict()[
                "_GasteigerCharge"
            ]
            charge_u = mappedmol.GetAtomWIthIdx(nb).GetPropsAsDict()["_GasteigerCharge"]
            if (
                mappedmol.GetAtomWithIdx(nbmap[nb][2]).GetSymbol()
                == mappedmol.GetAtomWIthIdx(nb).GetSymbol()
            ):  # No noticeable difference in electronegativity
                charge_u = 0
            elif charge_u < charge_m:
                charge_u = -1
            elif charge_u > charge_m:
                charge_u = 1
            else:
                raise Exception("Cannot determine charge of unmapped atom")
            fragmol = Chem.MolFromSmarts(candidate_strict)
            fragmol.UpdatePropertyCache(strict=False)

    except Exception as e:
        print(e)
        finallist += [candidatelist]
        oldfinallist += [oldcandidatelist]


def storeatommap(
    mappedsmiles, specid=0, idx=0, atommap={}, neighbormap={}, inverse=False
):
    if isinstance(mappedsmiles, str):
        mappedmol = molfromsmiles(mappedsmiles)
    else:
        mappedmol = mappedsmiles
    for atom in mappedmol.GetAtoms():
        if atom.HasProp("molAtomMapNumber"):
            mnum = atom.GetAtomMapNum()
            if not inverse:
                atommap[mnum] = (specid, idx, atom.GetIdx())
            else:
                atommap[atom.GetIdx()] = (
                    specid,
                    idx,
                    mnum,
                )  # NOTE: if aggregating across multiple species, inverse should be False otherwise information is rewritten
            #             breakpoint()
            neighbors = set(
                [
                    nb.GetIdx()
                    for nb in atom.GetNeighbors()
                    if not nb.HasProp("molAtomMapNumber")
                ]
            )
            if neighbors:
                if not inverse:
                    if neighbormap:
                        startidx = (
                            max([int(key.split("n")[1]) for key in neighbormap]) + 1
                        )
                    else:
                        startidx = 0
                    for i, nb in enumerate(neighbors, start=startidx):
                        neighbormap["n" + str(i)] = (specid, idx, nb, mnum)
                else:
                    for nb in neighbors:
                        neighbormap[nb] = (specid, idx, atom.GetIdx())
    return atommap, neighbormap


# def calcreliability()


def gen_template_fragment(
    fragidx,
    mappedmol,
    specificity="loose",
    atomsymbols=OrderedDict({}),
    funcgroupmapnum=set(),
    separateAllHs=False,
):
    #     breakpoint()
    if not atomsymbols:
        #         breakpoint()
        mappedmol = Chem.RemoveHs(mappedmol)  # Hydrogens become implicit
        clear_isotope(mappedmol)  # Redefining isotopes
        hatomsidx = {
            atom.GetIdx() for atom in mappedmol.GetAtoms() if atom.GetSymbol() == "H"
        }
        if fragidx.issubset(hatomsidx):
            hatomsremaining = {}
        else:
            hatomsremaining = fragidx.intersection(hatomsidx)
        if hatomsremaining:  # Likely isotopes
            affectedatomidx = {
                nb.GetIdx()
                for hatom in hatomsremaining
                for nb in mappedmol.GetAtomWithIdx(hatom).GetNeighbors()
            }
        else:
            affectedatomidx = {}
        fragidx = {
            idx
            for idx in fragidx
            if idx < len(mappedmol.GetAtoms())
            if idx not in hatomsremaining
        }
        for atom in Chem.AddHs(
            mappedmol
        ).GetAtoms():  # Functional group atoms described in detail, everything else has hydrogens stripped
            #             breakpoint()
            numHs_ = None
            degree_ = None
            if atom.GetSymbol() == "H" and atom.GetIdx() not in hatomsidx:
                continue
            refatom = mappedmol.GetAtomWithIdx(atom.GetIdx())
            if affectedatomidx and atom.GetIdx() in affectedatomidx:
                Hslist = [nb.GetSymbol() for nb in atom.GetNeighbors()]
                numHs_ = Hslist.count("H")
                degree_ = atom.GetDegree() - numHs_
            if specificity == "loose":
                if (
                    funcgroupmapnum
                    and atom.HasProp("molAtomMapNumber")
                    and atom.GetAtomMapNum() in funcgroupmapnum
                ):
                    #                 if refatom.GetDegree()==1: #Terminal atom
                    symbol = get_strict_smarts_for_atom(
                        refatom,
                        numHs_=numHs_,
                        degree_=degree_,
                        separateAllHs=separateAllHs,
                    )
                else:
                    symbol = atom.GetSmarts(
                        isomericSmiles=False
                    )  # General SMARTS to yield a generalizable template
            else:
                symbol = get_strict_smarts_for_atom(
                    refatom, numHs_=numHs_, degree_=degree_, separateAllHs=separateAllHs
                )
            atomsymbols.update({atom.GetIdx(): symbol})
    if specificity == "loose":
        frag = Chem.MolFragmentToSmiles(
            mappedmol, fragidx, atomSymbols=list(atomsymbols.values())
        )
    else:  # Closer to rdchiral template
        frag = Chem.MolFragmentToSmiles(
            mappedmol,
            fragidx,
            atomSymbols=list(atomsymbols.values()),
            allHsExplicit=True,
            allBondsExplicit=True,
        )
    return frag, atomsymbols, fragidx


def get_strict_smarts_for_atom(
    atom, numHs_=None, degree_=None, charge_=None, separateAllHs=False
):
    """
    For an RDkit atom object, generate a SMARTS pattern that
    matches the atom as strictly as possible, taken from rdChiral
    """
    USE_STEREOCHEMISTRY = True
    symbol = atom.GetSmarts()
    #     if atom.GetSymbol() == 'H':
    #         symbol = '[#1]'
    if numHs_ is not None:
        numHs = numHs_
    else:
        numHs = atom.GetTotalNumHs()
    if degree_ is not None:
        degree = degree_
    else:
        degree = atom.GetDegree()
    if charge_ is not None:
        charge = charge_
    else:
        charge = atom.GetFormalCharge()

    if "[" not in symbol:
        symbol = "[" + symbol + "]"

    # Explicit stereochemistry - *before* H
    if USE_STEREOCHEMISTRY:
        if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            if "@" not in symbol:
                # Be explicit when there is a tetrahedral chiral tag
                if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    tag = "@"
                elif atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                    tag = "@@"
                if ":" in symbol:
                    symbol = symbol.replace(":", ";{}:".format(tag))
                else:
                    symbol = symbol.replace("]", ";{}]".format(tag))

    if "H" not in symbol:
        H_symbol = "H{}".format(numHs)
        # Explicit number of hydrogens: include "H0" when no hydrogens present
        if ":" in symbol:  # stick H0 before label
            symbol = symbol.replace(":", ";{}:".format(H_symbol))
        else:
            symbol = symbol.replace("]", ";{}]".format(H_symbol))
    elif separateAllHs:
        if "H" in symbol and symbol != "[H]":
            smts = symbol.split("[")[1].split("]")[0].split(":")[0]
            parentatom = smts.split("H")[0]
            if parentatom:
                symbol = symbol.replace(smts, parentatom)
                H_symbol = "H{}".format(numHs)
                if ":" in symbol:  # stick H0 before label
                    symbol = symbol.replace(":", ";{}:".format(H_symbol))
                else:
                    symbol = symbol.replace("]", ";{}]".format(H_symbol))
    # Explicit degree
    if ":" in symbol:
        symbol = symbol.replace(":", ";D{}:".format(degree))
    else:
        symbol = symbol.replace("]", ";D{}]".format(degree))

    # Explicit formal charge
    if "+" not in symbol and "-" not in symbol:
        charge_symbol = "+" if (charge >= 0) else "-"
        charge_symbol += "{}".format(abs(charge))
        if ":" in symbol:
            symbol = symbol.replace(":", ";{}:".format(charge_symbol))
        else:
            symbol = symbol.replace("]", ";{}]".format(charge_symbol))

    return symbol


def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


# uspto_raw=pd.read_csv('/home/aa2133/Impurity-Project/USPTO/uspto_reactions_smiles_full.csv',index_col=0)
# rownumber=9400
# rxnsmiles0,Rdata,Pdata, balrxnsmiles, msg, LHSids, RHSids, hcrct, hcprod, LHSdata, RHSdata,Rgtdata,Solvdata, mappedrxn, conf, msg1=balance_rxn_uspto_row(uspto_raw.iloc[rownumber],hc_prod={},ignoreH=True)
# rxnsmiles=balrxnsmiles
# rxn=rdChemReactions.ReactionFromSmarts(mappedrxn)
# rxnsmiles02,balrxnsmiles2,msg0,LHSids2,RHSids2,hcrct2,hcprod2,LHSdata2,RHSdata2,Rgtdata2,Solvdata2,mappedrxn2,conf2,msg1_0=cleavagebal(LHSdata,RHSdata,rxnsmiles,rxn,msg=msg,Rdata=Rdata,Rgtdata=Rgtdata,Solvdata=Solvdata,resolvecharge=False) #Rdata=Rdata
# rxnsmiles02,balrxnsmiles2,msg0,LHSids2,RHSids2,hcrct2,hcprod2,LHSdata2,RHSdata2,_,_,mappedrxn2,conf2,msg1_0=resolvechargebal(LHSdata2,RHSdata2,Rdata,Rdata=Rdata,msg=msg0,Rgtdata=Rgtdata2,Solvdata=Solvdata2,ignoreH=True) #Rdata=Rdata
# balance_rxn_uspto_row(uspto_raw.iloc[1331])
# balance_rxn_dp("CCC(N)(CC)C(=O)OC.CCCOc1nc(C(=O)N[C@H](CO)CC(C)C)cnc1N1CCCC1.O=C(O)c1cnc(NC2CCCCC2)c(OCC(F)(F)F)n1>>CCC(CC)(NC(=O)c1cnc(NC2CCCCC2)c(OCC(F)(F)F)n1)C(=O)OC")
# demol=molfromsmiles("CC(=O)Oc1cc(OCc2ccccc2)ccc1OCc1ccccc1")
# gen_template_fragment({20,25,24},demol,specificity='strict',separateAllHs=False)
# rxnsmiles0_1,Rdata1,Pdata1, balrxnsmiles1, msg_1, LHSids1, RHSids1, hcrct1, hcprod1, LHSdata1, RHSdata1,Rgtdata1,Solvdata1, mappedrxn1, conf1, msg1_1=balance_rxn_uspto_row(uspto_raw.iloc[rownumber],ignoreH=True)
