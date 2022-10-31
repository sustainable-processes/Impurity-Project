from AnalgCompds import getCarrierFrags0
from MainFunctions import getfragments
from MapRxns import checkrxnrow, update_matches, updaterxns_
from RxnCenter import getrxncenterrow, validrxncenterrow
from helpCompound import hc_Dict
from typing import Dict, List
from rdkit import Chem  # Importing RDKit
from rdkit.Chem import rdChemReactions  # Reaction processing
from AnalgRxns import getspecdat_rxn
from BalanceRxns import balancerxn, maprxn
import pandas as pd
from GenTempl import gen_template_row


def balance_rxn_dp(rxn: str, **kwargs):
    """_summary_

    Args:
        rxn (str): _description_

    Returns:
        _type_: _description_
    """
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
        "first": False,
    }  # Help reactant dictionary
    IP = {**IP, **kwargs}
    Rdata, Pdata, Rgtdata, Solvdata = getspecdat_rxn(
        rxn, reagents=IP["reagents"], solvents=IP["solvents"]
    )
    IP["addedspecies"] = [i for i in Rdata]
    rxnsmiles0 = ">>".join(
        [
            getfragments([Rdata[r]["smiles"] for r in Rdata], smiles=True),
            getfragments([Pdata[p]["smiles"] for p in Pdata], smiles=True),
        ]
    )
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
    mappedrxn_ = maprxn([balrxnsmiles])[0]
    mappedrxn = mappedrxn_.get("mapped_rxn")
    conf = mappedrxn_.get("confidence")
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
        if "Unmapped" in msg1 or "unmapped" in msg1:
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
            ) = updaterxns_(
                pd.DataFrame(
                    [
                        {
                            "Rdata": Rdata,
                            "LHSdata": LHSdata,
                            "RHSdata": RHSdata,
                            "msg": msg,
                            "hcprod": hcprod,
                            "hcrct": hcrct,
                            "Rgtdata": Rgtdata,
                            "Solvdata": Solvdata,
                        }
                    ]
                ).iloc[0]
            )
    else:
        msg1 = "Mapping error"
    return (
        rxnsmiles0,
        balrxnsmiles,
        msg,
        LHSids,
        RHSids,
        hcrct,
        hcprod,
        LHSdata,
        RHSdata,
        mappedrxn,
        conf,
        msg1,
    )


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


def comparetemplate(reftemplate: str, template: str):
    reftemplrxn = rdChemReactions.ReactionFromSmarts(template)
    templrxn = rdChemReactions.ReactionFromSmarts(template)
