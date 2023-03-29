import functools
from collections import Counter
from typing import List, Union

import ipywidgets as widgets
import pandas as pd
from IPython.display import Markdown, clear_output, display
from ipywidgets import Layout
from rdkit import Chem  # Importing RDKit
from rdkit.Chem import Draw, rdChemReactions  # Reaction processing
from rdkit.Chem.Draw import IPythonConsole

from AnalgCompds import getCarrierFrags0
from FunctionsDB import getmixturefrags
from MainFunctions import (
    drawReaction,
    highlightsubstruct,
    molfromsmiles,
)
from MapRxns import maprxn
from Standalone import balance_rxn_dp, rxn_center_dp




# (
#     rxnsmiles0,
#     balrxnsmiles,
#     msg,
#     LHSids,
#     RHSids,
#     hcrct,
#     hcprod,
#     LHSdata,
#     RHSdata,
#     mappedrxn,
#     conf,
#     msg1,
# ) = balance_rxn_dp(
#     "Brc1csc2ccccc12.OB(O)c1ccc2c(c1)c1ccccc1n2-c1ccccc1>>c1ccc(-n2c3ccccc3c3cc(-c4csc5ccccc45)ccc32)cc1"
# )
# (
#     specmap,
#     rnbmap,
#     rxncentermapnum,
#     LHSdata,
#     msg,
#     outfrag,
#     outfg,
#     outneighbor,
#     unusedanalogue,
# ) = rxn_center_dp(mappedrxn, LHSdata, RHSdata)


def visfragment1(smiles: str, pattlist: List):
    return display(highlightsubstruct(smiles, pattlist=pattlist))


def verifyindex(
    df: pd.DataFrame, idxcol: Union[str, List[str]] = ["ReactionID", "Instance"]
):
    if isinstance(idxcol, list) and len(idxcol) == 1:
        idxcol = idxcol[0]
    if isinstance(idxcol, str):
        idx = df.index.name
        if idx and idx != idxcol:
            df.reset_index(inplace=True)
            df.set_index(idxcol, inplace=True)
        elif not idx:
            df.set_index(idxcol, inplace=True)
    else:
        idx = df.index.names
        if idx and idx != idxcol:
            df.reset_index(inplace=True)
            df.set_index(idxcol, inplace=True)
        elif not idx:
            df.set_index(idxcol, inplace=True)
    return df


def vismaster(**kwargs):
    # Setting up buttons
    dataminingbutton = widgets.Button(
        description="Data Mining",
        button_style="success",
        layout=widgets.Layout(width="20%", height="auto"),
    )
    dataprocessingbutton = widgets.Button(
        description="Data Processing",
        button_style="success",
        layout=widgets.Layout(width="20%", height="auto"),
    )
    impuritypredictionbutton = widgets.Button(
        description="Impurity Prediction",
        button_style="success",
        layout=widgets.Layout(width="20%", height="auto"),
    )
    impurityrankingbutton = widgets.Button(
        description="Impurity Ranking",
        button_style="success",
        layout=widgets.Layout(width="20%", height="auto"),
    )
    allbutton = widgets.Button(
        description="All",
        button_style="success",
        layout=widgets.Layout(width="25%", height="auto"),
    )
    uibutton = widgets.HBox(
        [
            allbutton,
            dataminingbutton,
            dataprocessingbutton,
            impuritypredictionbutton,
            impurityrankingbutton,
        ],
        layout=widgets.Layout(border="5px solid green"),
    )
    display(Markdown("#### <center>Select a stage to visualize</center>"))
    display(uibutton)
    outmaster = widgets.Output()

    def on_button_clicked(b, stages: List = []):
        with outmaster:
            clear_output()
            tracemaster(stages=stages, **kwargs)

    dataminingbutton.on_click(
        functools.partial(on_button_clicked, stages=["Data Mining"])
    )
    dataprocessingbutton.on_click(
        functools.partial(on_button_clicked, stages=["Data Processing"])
    )
    impuritypredictionbutton.on_click(
        functools.partial(on_button_clicked, stages=["Impurity Prediction"])
    )
    impurityrankingbutton.on_click(
        functools.partial(on_button_clicked, stages=["Impurity Ranking"])
    )
    allbutton.on_click(
        functools.partial(
            on_button_clicked,
            stages=[
                "Data Mining",
                "Data Processing",
                "Impurity Prediction",
                "Impurity Ranking",
            ],
        )
    )
    display(outmaster)


def tracemaster(
    stages=["Data Mining", "Data Processing", "Impurity Prediction"], **kwargs
):
    global displaywidget, optionschanged
    displaywidget = False
    optionschanged = False
    for stage in stages:
        if stage == "Data Mining":
            display(Markdown("<h1><center><strong>Data Mining</strong></center></h1>"))
            if "fragdbsource" in kwargs and isinstance(
                kwargs["fragdbsource"], pd.DataFrame
            ):
                fragdb = kwargs["fragdbsource"]
                fragdb = verifyindex(fragdb, idxcol=["FragmentSmiles", "SubstanceID"])
            else:
                fragdb = None
            for i, iq in enumerate(
                ["inputquery_analg_updated", "inputquery_analg", "inputquery"]
            ):
                if iq in kwargs and kwargs[iq] is not None:
                    inputquery = kwargs[iq]
                    userinput = inputquery["smiles"]
                    display(
                        Markdown(
                            f"<h2><center><strong>1. Query Reaction</strong></center></h2>"
                        )
                    )
                    display(
                        drawReaction(
                            rdChemReactions.ReactionFromSmarts(
                                userinput, useSmiles=True
                            )
                        )
                    )
                    display(Markdown(f"`Query Reaction SMILES: {userinput}`"))
                    queryspecies = list(inputquery["species"].keys())
                    queryspec = widgets.Dropdown(
                        options=[""] + queryspecies,
                        style={"description_width": "initial"},
                        description="Query Species",
                        continuous_update=False,
                    )
                    queryfrag = widgets.Dropdown(
                        description="Fragment", options=[""], continuous_update=False
                    )
                    queryexpand = widgets.IntSlider(
                        description="Expand",
                        min=0,
                        max=10,
                        step=1,
                        value=1,
                        continuous_update=False,
                    )
                    queryresformat = widgets.Dropdown(
                        options=["smiles", "smarts"],
                        style={"description_width": "initial"},
                        description="Result Format",
                        value="smiles",
                    )

                    def on_update_queryspec_widget(*args):
                        if queryspec.value:
                            if (
                                queryresformat.value == "smiles"
                                and queryexpand.value == 1
                            ):
                                queryfrags = list(
                                    inputquery["species"][queryspec.value].keys()
                                )
                            else:
                                if "." in queryspec.value:  # Mixture
                                    queryfrags = getmixturefrags(
                                        queryspec.value,
                                        expand=queryexpand.value,
                                        resFormat=queryresformat.value,
                                    )
                                else:
                                    queryfrags = getCarrierFrags0(
                                        queryspec.value,
                                        expand=queryexpand.value,
                                        resFormat=queryresformat.value,
                                    )
                            queryfrag.options = [""] + list(Counter(queryfrags).keys())
                        else:
                            queryfrag.options = [""]

                    def on_update_queryexpand_widget(*args):
                        if queryspec.value:
                            if (
                                queryexpand.value != 1
                                or queryresformat.value != "smiles"
                            ):
                                if "." in queryspec.value:  # Mixture
                                    queryfrags = getmixturefrags(
                                        queryspec.value,
                                        expand=queryexpand.value,
                                        resFormat=queryresformat.value,
                                    )
                                else:
                                    queryfrags = getCarrierFrags0(
                                        queryspec.value,
                                        expand=queryexpand.value,
                                        resFormat=queryresformat.value,
                                    )
                            else:
                                queryfrags = list(
                                    inputquery["species"][queryspec.value].keys()
                                )
                            queryfrag.options = [""] + list(Counter(queryfrags).keys())
                        else:
                            queryfrag.options = [""]

                    def on_update_queryresformat_widget(*args):
                        if queryspec.value:
                            if (
                                queryresformat.value != "smiles"
                                or queryexpand.value != 1
                            ):
                                if "." in queryspec.value:  # Mixture
                                    queryfrags = getmixturefrags(
                                        queryspec.value,
                                        expand=queryexpand.value,
                                        resFormat=queryresformat.value,
                                    )
                                else:
                                    queryfrags = getCarrierFrags0(
                                        queryspec.value,
                                        expand=queryexpand.value,
                                        resFormat=queryresformat.value,
                                    )
                            else:
                                queryfrags = list(
                                    inputquery["species"][queryspec.value].keys()
                                )
                            queryfrag.options = [""] + list(Counter(queryfrags).keys())
                        else:
                            queryfrag.options = [""]

                    def visfragment_(queryspec: str, frag: str):
                        if frag:
                            img, count = highlightsubstruct(
                                queryspec, [frag], returncount=True
                            )
                            display(img)
                            display(Markdown(f"Species SMILES: {queryspec}"))
                            display(Markdown(f"Fragment SMILES: {frag}"))
                            display(Markdown(f"Fragment count:{count[0]}"))

                            if i < 2:
                                specdata = inputquery["species"][queryspec]
                                if frag in specdata:
                                    analoguepool = specdata[frag]["analoguepool"]
                                    combinedui, fragui, outfrag = visfragment(
                                        analoguepool=analoguepool,
                                        analoguefrag=frag,
                                        displayres=False,
                                        **kwargs,
                                    )
                                    display(
                                        Markdown(
                                            "<h2><center><strong>3. Fragment Visualization (Analogue Species)</strong></center></h2>"
                                        )
                                    )
                                    display(
                                        Markdown(
                                            f"The chosen carrier fragment is {frag}"
                                        )
                                    )
                                    display(
                                        Markdown(
                                            "<h3><center><strong>Analogue Species</strong></center></h3>"
                                        )
                                    )
                                    display(combinedui)
                                    display(
                                        Markdown(
                                            "<h3><center><strong>Fragments</strong></center></h3>"
                                        )
                                    )
                                    display(fragui)
                                    display(outfrag)
                            else:
                                display(
                                    Markdown(
                                        f"<h3><center><strong>No analogue species information detected</strong></center></h3>"
                                    )
                                )

                    queryspec.observe(on_update_queryspec_widget, "value")
                    queryexpand.observe(on_update_queryexpand_widget, "value")
                    queryresformat.observe(on_update_queryresformat_widget, "value")
                    queryui_ = widgets.HBox([queryfrag, queryexpand, queryresformat])
                    queryui = widgets.VBox([queryspec, queryui_])
                    outmain = widgets.interactive_output(
                        visfragment_,
                        {
                            "queryspec": queryspec,
                            "frag": queryfrag,
                        },
                    )
                    display(
                        Markdown(
                            "<h2><center><strong>2. Fragment Visualization (Query Species)</strong></center></h2>"
                        )
                    )
                    display(queryui)

                    display(outmain)

                    break
        if stage == "Data Processing":
            displaywidget = False
            reactionid_dp = widgets.SelectionSlider(
                options=[""], description="Reaction ID", continuous_update=False
            )
            instance_dp = widgets.SelectionSlider(
                description="Instance", options=[0], continuous_update=False
            )
            reactionid_dptext = widgets.Text(description="Reaction ID", value="")
            customrxnsmiles_dp = widgets.Text(
                description="Custom Reaction SMILES",
                style={"description_width": "initial"},
            )
            # analgrawcheck = widgets.Checkbox(
            #     value=True,
            #     description="Restrict available reactions to analogue reactions from step 4",
            #     disabled=False,
            #     indent=False,
            # )
            # analgcheck = widgets.Checkbox(
            #     value=False,
            #     description="Restrict available reactions to this step",
            #     disabled=False,
            #     indent=False,
            # )
            # analgbalcheck = widgets.Checkbox(
            #     value=False,
            #     description="Restrict available reactions to this step",
            #     disabled=False,
            #     indent=False,
            # )
            # analgmapcheck = widgets.Checkbox(
            #     value=False,
            #     description="Restrict available reactions to this step",
            #     disabled=False,
            #     indent=False,
            # )
            # analgcentcheck = widgets.Checkbox(
            #     value=False,
            #     description="Restrict available reactions to this step",
            #     disabled=False,
            #     indent=False,
            # )

            # def on_update_analgcheckraw(*args):
            #     if analgrawcheck.value and kwargs["analoguerxns"] is not None:
            #         (
            #             reactionid_dp,
            #             reactionid_dptext,
            #             instance_dp,
            #             customrxnsmiles_dp,
            #         ) = updatereactionwidget(
            #             reactionid_dp, instance_dp, kwargs["analoguerxns"]
            #         )
            #         analgcheck.value = False
            #         analgbalcheck.value = False
            #         analgmapcheck.value = False
            #         analgcentcheck.value = False

            # def on_update_analgcheck(*args):
            #     if analgcheck.value and kwargs["analoguerxnsbal"] is not None:
            #         (
            #             reactionid_dp,
            #             reactionid_dptext,
            #             instance_dp,
            #             customrxnsmiles_dp,
            #         ) = updatereactionwidget(
            #             reactionid_dp, instance_dp, kwargs["analoguerxnsbal"]
            #         )
            #         analgrawcheck.value = False
            #         analgbalcheck.value = False
            #         analgmapcheck.value = False
            #         analgcentcheck.value = False
            #     elif all(
            #         [
            #             not i.value
            #             for i in [analgbalcheck, analgmapcheck, analgcentcheck]
            #         ]
            #     ):
            #         analgrawcheck.value = True

            # def on_update_analgbalcheck(*args):
            #     if analgbalcheck.value and kwargs["analoguerxnsmapped"] is not None:
            #         (
            #             reactionid_dp,
            #             reactionid_dptext,
            #             instance_dp,
            #             customrxnsmiles_dp,
            #         ) = updatereactionwidget(
            #             reactionid_dp, instance_dp, kwargs["analoguerxnsmapped"]
            #         )
            #         analgrawcheck.value = False
            #         analgcheck.value = False
            #         analgbalcheck.value = False
            #         analgcentcheck.value = False
            #     elif all(
            #         [not i.value for i in [analgcheck, analgmapcheck, analgcentcheck]]
            #     ):
            #         analgrawcheck.value = True

            # def on_update_analgmapcheck(*args):
            #     if analgmapcheck.value and kwargs["analoguerxnsassigned"] is not None:
            #         (
            #             reactionid_dp,
            #             reactionid_dptext,
            #             instance_dp,
            #             customrxnsmiles_dp,
            #         ) = updatereactionwidget(
            #             reactionid_dp, instance_dp, kwargs["analoguerxnsassigned"]
            #         )
            #         analgrawcheck.value = False
            #         analgcheck.value = False
            #         analgbalcheck.value = False
            #         analgcentcheck.value = False
            #     elif all(
            #         [not i.value for i in [analgcheck, analgbalcheck, analgcentcheck]]
            #     ):
            #         analgrawcheck.value = True

            # def on_update_analgcentcheck(*args):
            #     if analgcentcheck.value and kwargs["analoguerxnsfinal"] is not None:
            #         (
            #             reactionid_dp,
            #             reactionid_dptext,
            #             instance_dp,
            #             customrxnsmiles_dp,
            #         ) = updatereactionwidget(
            #             reactionid_dp, instance_dp, kwargs["analoguerxnsfinal"]
            #         )
            #         analgrawcheck.value = False
            #         analgcheck.value = False
            #         analgbalcheck.value = False
            #         analgmapcheck.value = False
            #     elif all(
            #         [not i.value for i in [analgcheck, analgbalcheck, analgmapcheck]]
            #     ):
            #         analgrawcheck.value = True

            # analgrawcheck.observe(on_update_analgcheckraw, "value")
            # analgcheck.observe(on_update_analgcheck, "value")
            # analgbalcheck.observe(on_update_analgbalcheck, "value")
            # analgmapcheck.observe(on_update_analgmapcheck, "value")
            # analgcentcheck.observe(on_update_analgcentcheck, "value")

            display(
                Markdown("<h1><center><strong>Data Processing</strong></center></h1>")
            )
            optionschanged = False
            for i, dfname in enumerate(
                [
                    "analoguerxns",
                    "analoguerxns_updated",
                    "analoguerxnsbal",
                    "analoguerxnsmapped",
                    "analoguerxnsparsed",
                    "analoguerxnsassigned",
                    "analoguerxnscent",
                    "analoguerxnsvalid",
                    "analoguerxnsfinal",
                ]
            ):  # All important inputs for data processing
                if dfname in kwargs:
                    df = kwargs[dfname]
                    if df is not None:
                        df = verifyindex(df)
                        if not optionschanged and i in [0, 2, 3, 5, 8]:
                            (
                                reactionid_dp,
                                reactionid_dptext,
                                instance_dp,
                                customrxnsmiles_dp,
                            ) = updatereactionwidget(
                                reactionid_dp,
                                instance_dp,
                                df,
                                reactionid_dptext,
                                customrxnsmiles_dp,
                            )
                            optionschanged = True
                            # if i == 0:
                            #     analgrawcheck.value = True
                            # elif i == 2:
                            #     analgcheck.value = True
                            # elif i == 3:
                            #     analgbalcheck.value = True
                            # elif i == 5:
                            #     analgmapcheck.value = True
                            # elif i == 8:
                            #     analgcentcheck.value = True
                else:
                    kwargs[dfname] = None

            def tracedprxns(
                reactionid: int,
                instance: int,
                customrxnsmiles: str,
                # analgrawcheck: bool,
                # analgcheck: bool,
                # analgbalcheck: bool,
                # analgmapcheck: bool,
                # analgcentcheck: bool,
            ):
                global displaywidget, success
                if kwargs["analoguerxns"] is not None:
                    display(
                        Markdown(
                            "<h2><center><strong>Analogue Reactions</strong></center></h2>"
                        )
                    )
                    display(
                        Markdown(
                            f"<center><strong>{len(kwargs['analoguerxns'])} reactions retrieved</strong></center>"
                        )
                    )
                reactionidui = widgets.VBox([reactionid_dp, reactionid_dptext])
                display(widgets.HBox([reactionidui, instance_dp, customrxnsmiles_dp]))
                displaywidget = True
                if (
                    customrxnsmiles
                ):  # User has defined a custom reaction SMILES (custom workflow is called)
                    try:
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
                            mappedrxn,
                            conf,
                            msg1,
                        ) = balance_rxn_dp(customrxnsmiles)
                        (
                            specmap,
                            rnbmap,
                            rxncentermapnum,
                            LHSdata,
                            msg,
                            outfrag,
                            outfg,
                            outneighbor,
                            unusedanalogue,
                        ) = rxn_center_dp(mappedrxn, LHSdata, RHSdata)
                    except Exception as e:
                        print(e)
                else:  # User can play around with existing reactions in supplied data frames
                    for i, dfname in enumerate(
                        [
                            "analoguerxnsfinal",
                            "analoguerxnsvalid",
                            "analoguerxnscent",
                            "analoguerxnsassigned",
                            "analoguerxnsparsed",
                            "analoguerxnsmapped",
                            "analoguerxnsbal",
                            "analoguerxns",
                        ]
                    ):
                        global success
                        success = False
                        if kwargs[dfname] is not None:
                            df = kwargs[dfname]
                            if (reactionid, instance) in df.index:
                                success = True
                                reactiondf = df.xs((reactionid, instance))
                                if i <= 7:
                                    if "rxnsmiles0" in reactiondf:
                                        display(
                                            Markdown(
                                                "<h2><center><strong>5. Analogue Reaction (Reaxys)</strong></center></h2>"
                                            )
                                        )
                                        rxninfo = {}
                                        if kwargs["analoguerxns_updated"] is not None:
                                            identifiers1 = [
                                                "Rdata",
                                                "Pdata",
                                                "Rgtdata",
                                                "Solvdata",
                                            ]
                                            rxninfo = (
                                                kwargs["analoguerxns_updated"]
                                                .xs((reactionid, instance))[
                                                    identifiers1
                                                ]
                                                .to_dict()
                                            )
                                        if kwargs["analoguerxns"] is not None:
                                            identifiers2 = [
                                                "MissingReactant",
                                                "MissingProduct",
                                                "MissingSolvent",
                                                "NameDict",
                                                "CatalystID",
                                                "MissingCatalyst",
                                                "Temperature",
                                                "Pressure",
                                                "ReactionTime",
                                                "YearPublished",
                                                "NumRefs",
                                                "NumSteps",
                                                "NumStages",
                                            ]
                                            rxninfo1 = (
                                                kwargs["analoguerxns"]
                                                .xs((reactionid, instance))[
                                                    identifiers2
                                                ]
                                                .to_dict()
                                            )
                                            rxninfo.update(rxninfo1)
                                        visreaction(reactiondf.rxnsmiles0, **rxninfo)
                                        if kwargs["analoguerxnsbal"] is not None:
                                            display(
                                                Markdown(
                                                    f"<center><strong>{len(kwargs['analoguerxnsbal'])} reactions remaining after filtering</strong></center>"
                                                )
                                            )
                                    if "balrxnsmiles" in reactiondf:
                                        display(
                                            Markdown(
                                                "<h2><center><strong>6. Balanced Reaction</strong></center></h2>"
                                            )
                                        )
                                        if reactiondf.balrxnsmiles != "Error":
                                            visreaction(reactiondf.balrxnsmiles)
                                        display(
                                            Markdown(
                                                f"<strong>Message: {reactiondf.msg}</center>"
                                            )
                                        )
                                        if kwargs["analoguerxnsmapped"] is not None:
                                            display(
                                                Markdown(
                                                    f"<center><strong>{len(kwargs['analoguerxnsmapped'])} reactions remaining after filtering</strong></center>"
                                                )
                                            )
                                    if "mapped_rxn" in reactiondf:
                                        display(
                                            Markdown(
                                                "<h2><center><strong>7. Mapped Reaction</strong></center></h2>"
                                            )
                                        )
                                        if reactiondf.mapped_rxn != "Error":
                                            visreaction(reactiondf.mapped_rxn)
                                            display(
                                                Markdown(
                                                    f"Confidence: {reactiondf.confidence}"
                                                )
                                            )
                                        else:
                                            display(Markdown("Error"))
                                if i <= 4:
                                    display(
                                        Markdown(
                                            f"<strong>Message: {reactiondf.msg1}</center>"
                                        )
                                    )
                                if i <= 3:
                                    if reactiondf.msg2 != reactiondf.msg1:
                                        display(
                                            Markdown(
                                                f"<strong>Message: {reactiondf.msg2}</center>"
                                            )
                                        )
                                    if kwargs["analoguerxnsassigned"] is not None:
                                        display(
                                            Markdown(
                                                f"<center><strong>{len(kwargs['analoguerxnsassigned'])} reactions remaining after filtering</strong></center>"
                                            )
                                        )
                                if i <= 2:
                                    if "rxncentermapnum" in reactiondf:
                                        display(
                                            Markdown(
                                                "<h2><center><strong>8. Reaction Center</strong></center></h2>"
                                            )
                                        )
                                        if reactiondf.rxncenter:
                                            display(
                                                Markdown(
                                                    f"<center><strong>Reaction center map numbers: {reactiondf.rxncentermapnum}</strong></center>"
                                                )
                                            )
                                        elif reactiondf.rnbmap:
                                            display(
                                                Markdown(
                                                    f"<center><strong>Reacting unmapped atoms: {reactiondf.rnbmap}</strong></center>"
                                                )
                                            )
                                        else:
                                            display(
                                                Markdown(
                                                    "<center><strong>Reaction does not have a reaction center</strong></center>"
                                                )
                                            )
                                IPythonConsole.drawOptions.setHighlightColour(
                                    (0.7, 1, 0.7)
                                )
                                IPythonConsole.drawOptions.minFontSize = 7
                                IPythonConsole.drawOptions.useBWAtomPalette()
                                IPythonConsole.drawOptions.legendFraction = 0.5
                                IPythonConsole.drawOptions.legendFontSize = 50
                                # IPythonConsole.drawOptions.legendFraction = 0.6
                                # IPythonConsole.drawOptions.useDefaultAtomPalette()

                                if i <= 1:
                                    colors = [(0.7, 1, 0.7), (1, 0.7, 0.7)]
                                    LHSdata = reactiondf.LHSdata
                                    drawsettings = {rctid: {} for rctid in LHSdata}
                                    for rctid in LHSdata:
                                        for ridx in LHSdata[rctid]["fragloc"]:
                                            highlightAtomList = [
                                                atm
                                                for frag in LHSdata[rctid]["fragloc"][
                                                    ridx
                                                ]
                                                for match in LHSdata[rctid]["fragloc"][
                                                    ridx
                                                ][frag]["corrmatches"]
                                                for atm in match
                                            ]
                                            if isinstance(ridx, tuple):
                                                rctmol = Chem.AddHs(
                                                    molfromsmiles(
                                                        LHSdata[rctid]["mappedsmiles"][
                                                            ridx[0]
                                                        ][ridx[1]]
                                                    )
                                                )
                                            else:
                                                rctmol = Chem.AddHs(
                                                    molfromsmiles(
                                                        LHSdata[rctid]["mappedsmiles"][
                                                            ridx
                                                        ]
                                                    )
                                                )
                                            highlightAtomColor = {
                                                atm: colors[0]
                                                for atm in highlightAtomList
                                            }
                                            legend = f"Species {rctid}, instance {ridx} \n Mapped SMILES: {LHSdata[rctid]['mappedsmiles'][ridx]}"
                                            if "reacfrag" in LHSdata[rctid] and LHSdata[
                                                rctid
                                            ]["reacfrag"].get(ridx):
                                                legend += f"\n Reacting carrier fragment(s): {', '.join(LHSdata[rctid]['reacfrag'][ridx].keys())}"
                                            drawsettings[rctid].update(
                                                {
                                                    ridx: [
                                                        rctmol,
                                                        highlightAtomList,
                                                        highlightAtomColor,
                                                        legend,
                                                        {},  # highlightBondColor
                                                    ]
                                                }
                                            )
                                    RHSdata = reactiondf.RHSdata
                                    drawsettingsp = {prodid: {} for prodid in RHSdata}
                                    for prodid in RHSdata:
                                        for i, mappedsmiles in enumerate(
                                            RHSdata[prodid]["mappedsmiles"]
                                        ):
                                            if isinstance(mappedsmiles, tuple):
                                                for j, mappedsmiles_ in enumerate(
                                                    mappedsmiles
                                                ):
                                                    pidx = (i, j)
                                                    prodmol = Chem.AddHs(
                                                        molfromsmiles(mappedsmiles_)
                                                    )
                                                    drawsettingsp[prodid].update(
                                                        {
                                                            pidx: [
                                                                prodmol,
                                                                [],
                                                                {},
                                                                f"Species {prodid}, instance {pidx} \n Mapped SMILES: {mappedsmiles_}",
                                                                {},
                                                            ]
                                                        }
                                                    )
                                            else:
                                                pidx = i
                                                prodmol = Chem.AddHs(
                                                    molfromsmiles(mappedsmiles)
                                                )
                                                drawsettingsp[prodid].update(
                                                    {
                                                        pidx: [
                                                            prodmol,
                                                            [],
                                                            {},
                                                            f"Species {prodid}, instance {pidx} \n Mapped SMILES: {mappedsmiles}",
                                                            {},
                                                        ]
                                                    }
                                                )

                                    RCs = [
                                        reactiondf.rxncentermapnum,
                                        set(reactiondf.rnbmap.keys()),
                                    ]
                                    mapdicts = [reactiondf.specmap, reactiondf.rnbmap]
                                    for RC, mapdict in zip(RCs, mapdicts):
                                        for changemapnum in RC:
                                            rctid = mapdict[changemapnum][0]
                                            ridx = mapdict[changemapnum][1]
                                            idxr = mapdict[changemapnum][2]
                                            prodid = mapdict[changemapnum][3]
                                            pidx = mapdict[changemapnum][4]
                                            idxp = mapdict[changemapnum][5]

                                            if (
                                                idxr not in drawsettings[rctid][ridx][1]
                                            ):  # Reaction center outside
                                                nbatoms = [
                                                    atom.GetIdx()
                                                    for atom in drawsettings[rctid][
                                                        ridx
                                                    ][0]
                                                    .GetAtomWithIdx(idxr)
                                                    .GetNeighbors()
                                                    if atom.GetIdx()
                                                    in drawsettings[rctid][ridx][1]
                                                ]
                                                if nbatoms:
                                                    drawsettings[rctid][ridx][4].update(
                                                        {
                                                            drawsettings[rctid][ridx][0]
                                                            .GetBondBetweenAtoms(
                                                                atomidx, idxr
                                                            )
                                                            .GetIdx(): colors[1]
                                                            for atomidx in nbatoms
                                                        }
                                                    )
                                                drawsettings[rctid][ridx][1].append(
                                                    idxr
                                                )
                                                drawsettings[rctid][ridx][2].update(
                                                    {idxr: colors[1]}
                                                )
                                            else:
                                                drawsettings[rctid][ridx][2][
                                                    idxr
                                                ] = colors[1]
                                            drawsettingsp[prodid][pidx][1].append(idxp)
                                            drawsettingsp[prodid][pidx][2].update(
                                                {idxp: colors[1]}
                                            )

                                    # print(drawsettings)
                                    # print(drawsettingsp)

                                    display(
                                        Markdown(
                                            "<h4><center><strong>Reactants</strong></center></h4>"
                                        )
                                    )
                                    display(
                                        Draw.MolsToGridImage(
                                            [
                                                drawsettings[rctid][ridx][0]
                                                for rctid in drawsettings
                                                for ridx in drawsettings[rctid]
                                            ],
                                            subImgSize=(500, 500),
                                            legends=[
                                                drawsettings[rctid][ridx][3]
                                                for rctid in drawsettings
                                                for ridx in drawsettings[rctid]
                                            ],
                                            highlightAtomLists=[
                                                drawsettings[rctid][ridx][1]
                                                for rctid in drawsettings
                                                for ridx in drawsettings[rctid]
                                            ],
                                            highlightAtomColors=[
                                                drawsettings[rctid][ridx][2]
                                                for rctid in drawsettings
                                                for ridx in drawsettings[rctid]
                                            ],
                                            highlightBondColors=[
                                                drawsettings[rctid][ridx][4]
                                                for rctid in drawsettings
                                                for ridx in drawsettings[rctid]
                                            ],
                                            useSVG=True,
                                        )
                                    )
                                    display(
                                        Markdown(
                                            f"<h4><center><strong>Products</strong></center></h4>"
                                        )
                                    )
                                    IPythonConsole.drawOptions.setHighlightColour(
                                        (1, 0.7, 0.7)
                                    )
                                    display(
                                        Draw.MolsToGridImage(
                                            [
                                                drawsettingsp[prodid][pidx][0]
                                                for prodid in drawsettingsp
                                                for pidx in drawsettingsp[prodid]
                                            ],
                                            subImgSize=(500, 500),
                                            legends=[
                                                drawsettingsp[prodid][pidx][3]
                                                for prodid in drawsettingsp
                                                for pidx in drawsettingsp[prodid]
                                            ],
                                            highlightAtomLists=[
                                                drawsettingsp[prodid][pidx][1]
                                                for prodid in drawsettingsp
                                                for pidx in drawsettingsp[prodid]
                                            ],
                                            highlightAtomColors=[
                                                drawsettingsp[prodid][pidx][2]
                                                for prodid in drawsettingsp
                                                for pidx in drawsettingsp[prodid]
                                            ],
                                            highlightBondColors=[
                                                drawsettingsp[prodid][pidx][4]
                                                for prodid in drawsettingsp
                                                for pidx in drawsettingsp[prodid]
                                            ],
                                            useSVG=True,
                                        )
                                    )

                                    display(
                                        Markdown(
                                            f"<center><strong>Message:{reactiondf.msg3}</strong></center>"
                                        )
                                    )
                                    if kwargs["analoguerxnsfinal"] is not None:
                                        display(
                                            Markdown(
                                                f"<center><strong>{len(kwargs['analoguerxnsfinal'])} reactions remaining after filtering</strong></center>"
                                            )
                                        )
                                break

            outdp = widgets.interactive_output(
                tracedprxns,
                {
                    "reactionid": reactionid_dp,
                    "instance": instance_dp,
                    "customrxnsmiles": customrxnsmiles_dp,
                    # "analgrawcheck": analgrawcheck,
                    # "analgcheck": analgcheck,
                    # "analgbalcheck": analgbalcheck,
                    # "analgmapcheck": analgmapcheck,
                    # "analgcentcheck": analgcentcheck,
                },
            )
            display(outdp)

        if stage == "Impurity Prediction":
            display(
                Markdown(
                    "<h1><center><strong>Impurity Prediction</strong></center></h1>"
                )
            )
            if not displaywidget:
                reactionid_ip = widgets.SelectionSlider(
                    options=[""], description="Reaction ID", continuous_update=False
                )
                instance_ip = widgets.SelectionSlider(
                    description="Instance", options=[0], continuous_update=False
                )
                reactionid_iptext = widgets.Text(description="Reaction ID", value="")
                customrxnsmiles_ip = widgets.Text(
                    description="Custom Reaction SMILES",
                    style={"description_width": "initial"},
                )
            else:
                reactionid_ip = reactionid_dp
                instance_ip = instance_dp
                reactionid_iptext = reactionid_dptext
                customrxnsmiles_ip = customrxnsmiles_dp

            for i, dfname in enumerate(
                [
                    "analoguerxnstempl",
                    "analoguerxnsimp",
                    "analoguerxnsimpfilt",
                    "impfinal",
                    "impfinalfilt",
                ]
            ):
                if dfname in kwargs:
                    df = kwargs[dfname]
                    if df is not None:
                        df = verifyindex(df)
                        if not optionschanged:
                            (
                                reactionid_ip,
                                reactionid_iptext,
                                instance_ip,
                                customrxnsmiles_ip,
                            ) = updatereactionwidget(
                                reactionid_ip,
                                instance_ip,
                                df,
                                reactionid_iptext,
                                customrxnsmiles_ip,
                            )
                            optionschanged = True
                else:
                    kwargs[dfname] = None

            def traceiprxns(reactionid, instance, customrxnsmiles):
                global displaywidget, success
                reactionidui = widgets.VBox([reactionid_ip, reactionid_iptext])
                display(widgets.HBox([reactionidui, instance_ip, customrxnsmiles_ip]))
                displaywidget = True

                for i, dfname in enumerate(
                    [
                        "impfinalfilt",
                        "impfinal",
                        "analoguerxnsimpfilt",
                        "analoguerxnsimp",
                        "analoguerxnstempl",
                    ]
                ):
                    success = False
                    if kwargs[dfname] is not None:
                        df = kwargs[dfname]
                        if (reactionid, instance) in df.index:
                            success = True
                            reactiondf = df.xs((reactionid, instance))
                            if i <= 4:
                                if "template" in reactiondf:
                                    if isinstance(reactiondf, pd.DataFrame):
                                        reactiondf_ = reactiondf.iloc[0]
                                    else:
                                        reactiondf_ = reactiondf
                                    display(
                                        Markdown(
                                            f"<h2><center><strong>9. Template</strong></center></h2>"
                                        )
                                    )
                                    display(
                                        drawReaction(
                                            rdChemReactions.ReactionFromSmarts(
                                                reactiondf_.template
                                            )
                                        )
                                    )
                                    display(
                                        Markdown(
                                            f"`Template SMARTS: {reactiondf_.template}`"
                                        )
                                    )
                                    display(
                                        Markdown(
                                            f"<strong>Message:{reactiondf_.msg4}</strong>"
                                        )
                                    )
                            if i <= 3:
                                if kwargs["analoguerxnsimp"] is not None and i != 3:
                                    reactiondf_ = kwargs["analoguerxnsimp"].xs(
                                        (reactionid, instance)
                                    )
                                else:
                                    reactiondf_ = reactiondf
                                if "impurityrxn" in reactiondf_:
                                    display(
                                        Markdown(
                                            f"<h2><center><strong>10. Template Application</strong></center></h2>"
                                        )
                                    )
                                    displaymsg = False
                                    if isinstance(reactiondf_, pd.DataFrame):
                                        for idx, reactiondf__ in reactiondf.iterrows():
                                            if (
                                                reactiondf__.impurityrxn
                                                and reactiondf__.impurityrxn != "Error"
                                            ):
                                                display(
                                                    drawReaction(
                                                        rdChemReactions.ReactionFromSmarts(
                                                            reactiondf__.impurityrxn,
                                                            useSmiles=True,
                                                        )
                                                    )
                                                )
                                                display(
                                                    Markdown(
                                                        f"`Impurity Reaction SMILES: {reactiondf__.impurityrxn}`"
                                                    )
                                                )
                                            if not displaymsg:
                                                display(
                                                    Markdown(
                                                        f"<strong>Message:{reactiondf__.msg5}</strong>"
                                                    )
                                                )
                                                displaymsg = True
                                    else:
                                        if (
                                            reactiondf_.impurityrxn
                                            and reactiondf_.impurityrxn != "Error"
                                        ):
                                            for impurityrxn in reactiondf_.impurityrxn:
                                                display(
                                                    drawReaction(
                                                        rdChemReactions.ReactionFromSmarts(
                                                            list(impurityrxn)[0],
                                                            useSmiles=True,
                                                        )
                                                    )
                                                )
                                                display(
                                                    Markdown(
                                                        f"`Impurity Reaction SMILES: {list(impurityrxn)[0]}`"
                                                    )
                                                )
                                        display(
                                            Markdown(
                                                f"<strong>Message:{reactiondf_.msg5}</strong>"
                                            )
                                        )
                                if kwargs["impfinal"] is not None:
                                    display(
                                        Markdown(
                                            f"<center><strong>{len(kwargs['impfinal'].index.unique())} reactions remaining after filtering</strong></center>"
                                        )
                                    )
                            if i <= 1:
                                if "msg6" in reactiondf:
                                    display(
                                        Markdown(
                                            f"<h2><center><strong>11. Impurity Cleaning</strong></center></h2>"
                                        )
                                    )
                                    if isinstance(reactiondf, pd.DataFrame):
                                        reactiondf = reactiondf.iloc[0]
                                    if isinstance(reactiondf.msg6, pd.DataFrame):
                                        msg6 = reactiondf.msg6.iloc[0]
                                    else:
                                        msg6 = reactiondf.msg6
                                    display(
                                        Markdown(f"<strong>Message:{msg6}</strong>")
                                    )
                                if kwargs["impfinalfilt"] is not None:
                                    display(
                                        Markdown(
                                            f"<center><strong>{len(kwargs['impfinalfilt'].index.unique())} reactions remaining after filtering</strong></center>"
                                        )
                                    )

                            break

            outip = widgets.interactive_output(
                traceiprxns,
                {
                    "reactionid": reactionid_ip,
                    "instance": instance_ip,
                    "customrxnsmiles": customrxnsmiles_ip,
                    # "analgrawcheck": analgrawcheck,
                    # "analgcheck": analgcheck,
                    # "analgbalcheck": analgbalcheck,
                    # "analgmapcheck": analgmapcheck,
                    # "analgcentcheck": analgcentcheck,
                },
            )

            display(outip)

        if stage == "Impurity Ranking":
            display(
                Markdown("<h1><center><strong>Impurity Ranking</strong></center></h1>")
            )
            for i, iq in enumerate(
                ["inputquery_analg_updated", "inputquery_analg", "inputquery"]
            ):
                if iq in kwargs and kwargs[iq] is not None:
                    inputquery = kwargs[iq]
                    userinput = inputquery["smiles"]
                    display(
                        Markdown(
                            f"<h2><center><strong>Query Reaction (User Input)</strong></center></h2>"
                        )
                    )
                    display(
                        drawReaction(
                            rdChemReactions.ReactionFromSmarts(
                                userinput, useSmiles=True
                            )
                        )
                    )
                    display(Markdown(f"`Query Reaction SMILES: {userinput}`"))
                    break
            for i, dfname in enumerate(
                [
                    "impfinal",
                    "impfinalfilt",
                    "impfinalfilt2",
                    "impfinalfilt3",
                    "impfinalfilt4",
                    "impfinalfilt5",
                ]
            ):
                if dfname in kwargs:
                    df = kwargs[dfname]
                    if df is not None:
                        df = verifyindex(df)
                else:
                    kwargs[dfname] = None
            impview = widgets.Button(
                description="Impurity View",
                button_style="success",
                layout=widgets.Layout(width="50%", height="auto"),
            )
            rxnview = widgets.Button(
                description="Reaction View",
                button_style="success",
                layout=widgets.Layout(width="50%", height="auto"),
            )
            iruibutton = widgets.HBox(
                [impview, rxnview], layout=widgets.Layout(border="5px solid green")
            )
            display(Markdown("## <center>Select a view</center>"))
            display(iruibutton)
            irmaster = widgets.Output()

            def on_view_clicked(b, view="Impurity"):
                with irmaster:
                    clear_output()
                    irview(view=view, **kwargs)

            def irview(view="Impurity", **kwargs):
                if kwargs["summary3"] is not None:
                    display(
                        Markdown(
                            "<h2><strong><center>Main Product (Query) Reaction</center></strong></h2>"
                        )
                    )
                    visreaction(kwargs["summary3"].iloc[0].rxn)
                    display(
                        Markdown(
                            f"Max relevance: {kwargs['summary3'].iloc[0]['''Max relevance_tfiltered''']}"
                        )
                    )
                    display(
                        Markdown(
                            f"Number of hits: {kwargs['summary3'].iloc[0].Hits_tfiltered}"
                        )
                    )
                    Trange = kwargs["summary3"].iloc[0].t_range
                    display(
                        Markdown(
                            f"Temperature range (&deg;C): {'-'.join([str(temp) for temp in Trange])}"
                        )
                    )
                    display(
                        Markdown(
                            "*Temperature range is based on 5th to 95th percentile of top 10 % relevant reactions*"
                        )
                    )
                    errorcodedict = {
                        "impfinalfilt2": "Invalid condition",
                        "impfinalfilt3": "Suspect self-reaction (without catalyst or reagent)",
                        "impfinalfilt4": "Invalid catalyst",
                        "impfinalfilt5": "Missing reactants/products",
                        "Frame2": "Missing temperature",
                        "Frame4": "Temperature outside indicated range",
                    }
                    if view == "Impurity":  # Impurity view
                        mainprods = kwargs["summary3"].iloc[0].products
                        impcountwidget = widgets.IntSlider(
                            description="Number of impurity reactions to display",
                            style={"description_width": "initial"},
                            layout=Layout(width="500px"),
                            min=1,
                            max=5,
                            step=1,
                            value=2,
                        )

                        def displayimpurities(impuritycount):
                            display(
                                Markdown(
                                    f"## <center><strong>Impurity Reactions</strong></center>"
                                )
                            )
                            display(impcountwidget)
                            i = 1
                            j = 1
                            while i <= impuritycount:
                                imprxn = kwargs["summary3"].iloc[j].rxn
                                Trangei = kwargs["summary3"].iloc[j].t_range
                                if min(Trangei) > max(Trange) or max(Trangei) < min(
                                    Trange
                                ):
                                    j += 1
                                    continue
                                impprods = kwargs["summary3"].iloc[j].products
                                if (
                                    "hc_Dict" in kwargs
                                    and kwargs["hc_Dict"] is not None
                                ):
                                    if i > 0 and any(
                                        [
                                            impprod in mainprods
                                            and impprod
                                            not in [
                                                kwargs["hc_Dict"][k]["smiles"]
                                                for k in kwargs["hc_Dict"]
                                            ]
                                            for impprod in impprods
                                        ]
                                    ):
                                        j += 1
                                        continue

                                display(
                                    Markdown(
                                        f"<h2><strong><center>{i}.</center></strong></h2>"
                                    )
                                )
                                visreaction(imprxn)
                                display(
                                    Markdown(
                                        f"Max relevance: {kwargs['summary3'].iloc[j]['''Max relevance_tfiltered''']}"
                                    )
                                )
                                display(
                                    Markdown(
                                        f"Number of hits: {kwargs['summary3'].iloc[j].Hits_tfiltered}"
                                    )
                                )
                                display(
                                    Markdown(
                                        f"Temperature range (&deg;C): {'-'.join([str(temp) for temp in Trangei])}"
                                    )
                                )
                                display(
                                    Markdown(
                                        "*Temperature range is based on 5th to 95th percentile of top 10 % relevant reactions*"
                                    )
                                )
                                reactionid_imp = widgets.SelectionSlider(
                                    options=[""],
                                    description="Reaction ID",
                                    continuous_update=False,
                                )
                                instance_imp = widgets.SelectionSlider(
                                    description="Instance",
                                    options=[0],
                                    continuous_update=False,
                                )
                                reactionid_imptext = widgets.Text(
                                    description="Reaction ID", value=""
                                )

                                customrxnsmiles_imp = widgets.Text(
                                    description="Custom Reaction SMILES",
                                    style={"description_width": "initial"},
                                )

                                (
                                    reactionid_imp,
                                    reactionid_imptext,
                                    instance_imp,
                                    customrxnsmiles_imp,
                                ) = updatereactionwidget(
                                    reactionid_imp,
                                    instance_imp,
                                    kwargs["summary3"].iloc[j].Frame4,
                                    reactionid_imptext,
                                    customrxnsmiles_imp,
                                )
                                display(
                                    Markdown("## <center>Analogue Reactions</center>")
                                )
                                impui = widgets.VBox(
                                    [reactionid_imp, reactionid_imptext]
                                )
                                display(
                                    widgets.HBox(
                                        [impui, instance_imp],
                                        layout=widgets.Layout(border="5px solid green"),
                                    )
                                )
                                analgframe = kwargs["summary3"].iloc[j].Frame4
                                i += 1
                                j += 1

                                def traceimprxns(reactionid, instance):
                                    if (reactionid, instance) in analgframe.index:
                                        reactiondf = impframe_old.xs(
                                            (reactionid, instance)
                                        )
                                        if isinstance(reactiondf, pd.DataFrame):
                                            reactiondf = reactiondf.iloc[0]
                                        display(
                                            drawReaction(
                                                rdChemReactions.ReactionFromSmarts(
                                                    reactiondf.mapped_rxn,
                                                    useSmiles=True,
                                                )
                                            )
                                        )
                                        display(
                                            Markdown(
                                                f"`Reaction SMILES: {reactiondf.mapped_rxn}`"
                                            )
                                        )
                                        display(
                                            Markdown(
                                                f"Relevance: {reactiondf.Relevance_morgan}"
                                            )
                                        )
                                        LHSdata = {
                                            specid: {
                                                "smiles": reactiondf.LHSdata[specid][
                                                    "smiles"
                                                ],
                                                "count": reactiondf.LHSdata[specid][
                                                    "count"
                                                ],
                                                "name": reactiondf.NameDict[specid],
                                            }
                                            for specid in reactiondf.LHSdata
                                        }
                                        display(Markdown(f"Reactant Data: {LHSdata}"))
                                        rgts = set(reactiondf.ReagentID) - set(
                                            reactiondf.LHS
                                        )
                                        if rgts:
                                            Rgtdata = {
                                                specid: {
                                                    "smiles": reactiondf.Rgtdata[
                                                        specid
                                                    ]["smiles"],
                                                    "count": reactiondf.Rgtdata[specid][
                                                        "count"
                                                    ],
                                                    "name": reactiondf.NameDict[specid],
                                                }
                                                for specid in rgts
                                            }
                                            display(
                                                Markdown(f"Reagent Data: {Rgtdata}")
                                            )

                                        cats = reactiondf.CatalystID2
                                        if cats:
                                            Catdata = {
                                                specid: {
                                                    "name": reactiondf.NameDict[specid]
                                                }
                                                for specid in cats
                                            }
                                            display(
                                                Markdown(f"Catalyst Data: {Catdata}")
                                            )
                                        solvs = reactiondf.SolventID
                                        if solvs:
                                            if "Solvdata" in reactiondf:
                                                Solvdata = {
                                                    specid: {
                                                        "smiles": reactiondf.Solvdata[
                                                            specid
                                                        ]["smiles"],
                                                        "count": reactiondf.Solvdata[
                                                            specid
                                                        ]["count"],
                                                        "name": reactiondf.NameDict[
                                                            specid
                                                        ],
                                                    }
                                                    for specid in solvs
                                                }
                                            elif "analoguerxns_updated" in kwargs:
                                                Solvdata = {
                                                    specid: {
                                                        "smiles": kwargs[
                                                            "analoguerxns_updated"
                                                        ]
                                                        .xs((reactionid, instance))
                                                        .Solvdata[specid]["smiles"],
                                                        "count": kwargs[
                                                            "analoguerxns_updated"
                                                        ]
                                                        .xs((reactionid, instance))
                                                        .Solvdata[specid]["count"],
                                                        "name": reactiondf.NameDict[
                                                            specid
                                                        ],
                                                    }
                                                    for specid in solvs
                                                }
                                            else:
                                                Solvdata = {
                                                    specid: {
                                                        "name": reactiondf.NameDict[
                                                            specid
                                                        ]
                                                    }
                                                    for specid in solvs
                                                }
                                            display(
                                                Markdown(f"Solvent Data: {Solvdata}")
                                            )
                                        RHSdata = {
                                            specid: {
                                                "smiles": reactiondf.RHSdata[specid][
                                                    "smiles"
                                                ],
                                                "count": reactiondf.RHSdata[specid][
                                                    "count"
                                                ],
                                            }
                                            for specid in reactiondf.RHSdata
                                        }
                                        for specid in RHSdata:
                                            if specid not in set(reactiondf.hcprod):
                                                RHSdata[specid].update(
                                                    {
                                                        "name": reactiondf.NameDict[
                                                            specid
                                                        ]
                                                    }
                                                )
                                        display(Markdown(f"Product Data: {RHSdata}"))
                                        if reactiondf.Temperature:
                                            display(
                                                Markdown(
                                                    f"Temperature (&deg;C): {reactiondf.Temperature}"
                                                )
                                            )
                                        if reactiondf.Pressure:
                                            display(
                                                Markdown(
                                                    f"Pressure (bar): {reactiondf.Pressure}"
                                                )
                                            )
                                        if reactiondf.ReactionTime:
                                            display(
                                                Markdown(
                                                    f"Reaction Time (hours): {reactiondf.ReactionTime}"
                                                )
                                            )
                                        display(
                                            Markdown(
                                                f"Number of References:{reactiondf.NumRefs}"
                                            )
                                        )
                                        display(
                                            Markdown(
                                                f"Number of Steps:{reactiondf.NumSteps}"
                                            )
                                        )
                                        display(
                                            Markdown(
                                                f"Number of Stages:{reactiondf.NumStages}"
                                            )
                                        )
                                        if reactiondf.YearPublished:
                                            display(
                                                Markdown(
                                                    f"Year Published:{reactiondf.YearPublished}"
                                                )
                                            )
                                        if reactiondf.ConditionNotes:
                                            display(
                                                Markdown(
                                                    f"Condition Notes:{reactiondf.ConditionNotes}"
                                                )
                                            )

                                outimpa = widgets.interactive_output(
                                    traceimprxns,
                                    {
                                        "reactionid": reactionid_imp,
                                        "instance": instance_imp,
                                    },
                                )
                                display(outimpa)

                                if kwargs["impfinalfilt"] is not None:
                                    impframe_old = (
                                        kwargs["impfinalfilt"]
                                        .loc[
                                            kwargs["impfinalfilt"].impurityrxn == imprxn
                                        ]
                                        .sort_values(
                                            by="Relevance_morgan", ascending=False
                                        )
                                    )
                                    rejrxns = impframe_old.loc[
                                        impframe_old.Relevance_morgan
                                        > kwargs["summary3"].iloc[i][
                                            """Max relevance_tfiltered"""
                                        ]
                                    ]
                                    if not rejrxns.empty:
                                        reactionid_rej = widgets.SelectionSlider(
                                            options=[""],
                                            description="Reaction ID",
                                            continuous_update=False,
                                        )
                                        instance_rej = widgets.SelectionSlider(
                                            description="Instance",
                                            options=[0],
                                            continuous_update=False,
                                        )
                                        reactionid_rejtext = widgets.Text(
                                            description="Reaction ID", value=""
                                        )
                                        customrxnsmiles_rej = widgets.Text(
                                            description="Custom Reaction SMILES",
                                            style={"description_width": "initial"},
                                        )

                                        (
                                            reactionid_rej,
                                            reactionid_rejtext,
                                            instance_rej,
                                            customrxnsmiles_rej,
                                        ) = updatereactionwidget(
                                            reactionid_rej,
                                            instance_rej,
                                            rejrxns,
                                            reactionid_rejtext,
                                            customrxnsmiles_rej,
                                        )

                                        def tracerejrxns(reactionid, instance):
                                            if (
                                                reactionid,
                                                instance,
                                            ) in rejrxns.index:
                                                reactiondf = rejrxns.xs(
                                                    (reactionid, instance)
                                                )
                                                if isinstance(reactiondf, pd.DataFrame):
                                                    reactiondf = reactiondf.iloc[0]
                                                display(
                                                    drawReaction(
                                                        rdChemReactions.ReactionFromSmarts(
                                                            reactiondf.mapped_rxn,
                                                            useSmiles=True,
                                                        )
                                                    )
                                                )
                                                display(
                                                    Markdown(
                                                        f"`Reaction SMILES: {reactiondf.mapped_rxn}`"
                                                    )
                                                )
                                                for dfname in errorcodedict:
                                                    if (
                                                        dfname in kwargs
                                                        and kwargs[dfname] is not None
                                                    ):
                                                        if (
                                                            reactionid,
                                                            instance,
                                                        ) not in kwargs[dfname].index:
                                                            errorcode = errorcodedict[
                                                                dfname
                                                            ]
                                                            display(
                                                                Markdown(
                                                                    f"<strong>Reaction removed due to: {errorcode}</strong>"
                                                                )
                                                            )
                                                            break
                                                    elif dfname == "Frame2":
                                                        if (
                                                            reactionid,
                                                            instance,
                                                        ) not in kwargs[
                                                            "summary3"
                                                        ].iloc[
                                                            i
                                                        ].Frame3.index:
                                                            errorcode = errorcodedict[
                                                                dfname
                                                            ]
                                                            display(
                                                                Markdown(
                                                                    f"Reaction removed due to: {errorcode}"
                                                                )
                                                            )
                                                            break
                                                    elif dfname == "Frame4":
                                                        if (
                                                            reactionid,
                                                            instance,
                                                        ) not in kwargs[
                                                            "summary3"
                                                        ].iloc[
                                                            i
                                                        ].Frame4.index:
                                                            errorcode = errorcodedict[
                                                                dfname
                                                            ]
                                                            display(
                                                                Markdown(
                                                                    f"Reaction removed due to: {errorcode}"
                                                                )
                                                            )
                                                            break
                                                display(
                                                    Markdown(
                                                        f"Relevance: {reactiondf.Relevance_morgan}"
                                                    )
                                                )
                                                LHSdata = {
                                                    specid: {
                                                        "smiles": reactiondf.LHSdata[
                                                            specid
                                                        ]["smiles"],
                                                        "count": reactiondf.LHSdata[
                                                            specid
                                                        ]["count"],
                                                        "name": reactiondf.NameDict[
                                                            specid
                                                        ],
                                                    }
                                                    for specid in reactiondf.LHSdata
                                                }
                                                display(
                                                    Markdown(
                                                        f"Reactant Data: {LHSdata}"
                                                    )
                                                )
                                                rgts = set(reactiondf.ReagentID) - set(
                                                    reactiondf.LHS
                                                )
                                                if rgts:
                                                    Rgtdata = {
                                                        specid: {
                                                            "smiles": reactiondf.Rgtdata[
                                                                specid
                                                            ][
                                                                "smiles"
                                                            ],
                                                            "count": reactiondf.Rgtdata[
                                                                specid
                                                            ]["count"],
                                                            "name": reactiondf.NameDict[
                                                                specid
                                                            ],
                                                        }
                                                        for specid in rgts
                                                    }
                                                    display(
                                                        Markdown(
                                                            f"Reagent Data: {Rgtdata}"
                                                        )
                                                    )

                                                cats = reactiondf.CatalystID2
                                                if cats:
                                                    Catdata = {
                                                        specid: {
                                                            "name": reactiondf.NameDict[
                                                                specid
                                                            ]
                                                        }
                                                        for specid in cats
                                                    }
                                                    display(
                                                        Markdown(
                                                            f"Catalyst Data: {Catdata}"
                                                        )
                                                    )
                                                missingcatalyst = (
                                                    reactiondf.MissingCatalyst
                                                )
                                                if missingcatalyst:
                                                    display(
                                                        Markdown(
                                                            f"Missing Catalyst: {missingcatalyst}"
                                                        )
                                                    )
                                                solvs = reactiondf.SolventID
                                                if solvs:
                                                    if "Solvdata" in reactiondf:
                                                        Solvdata = {
                                                            specid: {
                                                                "smiles": reactiondf.Solvdata[
                                                                    specid
                                                                ][
                                                                    "smiles"
                                                                ],
                                                                "count": reactiondf.Solvdata[
                                                                    specid
                                                                ][
                                                                    "count"
                                                                ],
                                                                "name": reactiondf.NameDict[
                                                                    specid
                                                                ],
                                                            }
                                                            for specid in solvs
                                                        }
                                                    elif (
                                                        "analoguerxns_updated" in kwargs
                                                    ):
                                                        Solvdata = {
                                                            specid: {
                                                                "smiles": kwargs[
                                                                    "analoguerxns_updated"
                                                                ]
                                                                .xs(
                                                                    (
                                                                        reactionid,
                                                                        instance,
                                                                    )
                                                                )
                                                                .Solvdata[specid][
                                                                    "smiles"
                                                                ],
                                                                "count": kwargs[
                                                                    "analoguerxns_updated"
                                                                ]
                                                                .xs(
                                                                    (
                                                                        reactionid,
                                                                        instance,
                                                                    )
                                                                )
                                                                .Solvdata[specid][
                                                                    "count"
                                                                ],
                                                                "name": reactiondf.NameDict[
                                                                    specid
                                                                ],
                                                            }
                                                            for specid in solvs
                                                        }
                                                    else:
                                                        Solvdata = {
                                                            specid: {
                                                                "name": reactiondf.NameDict[
                                                                    specid
                                                                ]
                                                            }
                                                            for specid in solvs
                                                        }
                                                    display(
                                                        Markdown(
                                                            f"Solvent Data: {Solvdata}"
                                                        )
                                                    )
                                                missingsolvent = (
                                                    reactiondf.MissingSolvent
                                                )
                                                if missingsolvent:
                                                    display(
                                                        Markdown(
                                                            f"Missing Solvent: {reactiondf.MissingSolvent}"
                                                        )
                                                    )
                                                RHSdata = {
                                                    specid: {
                                                        "smiles": reactiondf.RHSdata[
                                                            specid
                                                        ]["smiles"],
                                                        "count": reactiondf.RHSdata[
                                                            specid
                                                        ]["count"],
                                                    }
                                                    for specid in reactiondf.RHSdata
                                                }
                                                for specid in RHSdata:
                                                    if specid not in set(
                                                        reactiondf.hcprod
                                                    ):
                                                        RHSdata[specid].update(
                                                            {
                                                                "name": reactiondf.NameDict[
                                                                    specid
                                                                ]
                                                            }
                                                        )
                                                display(
                                                    Markdown(f"Product Data: {RHSdata}")
                                                )
                                                if reactiondf.Temperature:
                                                    display(
                                                        Markdown(
                                                            f"Temperature (&deg;C): {reactiondf.Temperature}"
                                                        )
                                                    )
                                                if reactiondf.Pressure:
                                                    display(
                                                        Markdown(
                                                            f"Pressure (bar): {reactiondf.Pressure}"
                                                        )
                                                    )
                                                if reactiondf.ReactionTime:
                                                    display(
                                                        Markdown(
                                                            f"Reaction Time (hours): {reactiondf.ReactionTime}"
                                                        )
                                                    )
                                                display(
                                                    Markdown(
                                                        f"Number of References:{reactiondf.NumRefs}"
                                                    )
                                                )
                                                display(
                                                    Markdown(
                                                        f"Number of Steps:{reactiondf.NumSteps}"
                                                    )
                                                )
                                                display(
                                                    Markdown(
                                                        f"Number of Stages:{reactiondf.NumStages}"
                                                    )
                                                )
                                                if reactiondf.YearPublished:
                                                    display(
                                                        Markdown(
                                                            f"Year Published:{reactiondf.YearPublished}"
                                                        )
                                                    )
                                                if reactiondf.ConditionNotes:
                                                    display(
                                                        Markdown(
                                                            f"Condition Notes:{reactiondf.ConditionNotes}"
                                                        )
                                                    )

                                        outrej = widgets.interactive_output(
                                            tracerejrxns,
                                            {
                                                "reactionid": reactionid_rej,
                                                "instance": instance_rej,
                                            },
                                        )

                                        display(
                                            Markdown(
                                                f"## <center>Rejected Analogue Reactions (>{kwargs['summary3'].iloc[i]['''Max relevance_tfiltered''']} relevance)</center>"
                                            )
                                        )
                                        rejui = widgets.VBox(
                                            [reactionid_rej, reactionid_rejtext]
                                        )
                                        display(
                                            widgets.HBox(
                                                [rejui, instance_rej],
                                                layout=widgets.Layout(
                                                    border="5px solid green"
                                                ),
                                            )
                                        )
                                        display(outrej)

                        outimp = widgets.interactive_output(
                            displayimpurities, {"impuritycount": impcountwidget}
                        )
                        display(outimp)

                    elif view == "Reaction":
                        if not displaywidget:
                            reactionid_ir = widgets.SelectionSlider(
                                options=[""],
                                description="Reaction ID",
                                continuous_update=False,
                            )
                            instance_ir = widgets.SelectionSlider(
                                description="Instance",
                                options=[0],
                                continuous_update=False,
                            )
                            reactionid_irtext = widgets.Text(
                                description="Reaction ID", value=""
                            )
                            customrxnsmiles_ir = widgets.Text(
                                description="Custom Reaction SMILES",
                                style={"description_width": "initial"},
                            )
                        else:
                            reactionid_ir = reactionid_dp
                            instance_ir = instance_dp
                            reactionid_irtext = reactionid_dptext
                            customrxnsmiles_ir = customrxnsmiles_dp
                        for i, dfname in enumerate(
                            [
                                "impfinal",
                                "impfinalfilt",
                                "impfinalfilt2",
                                "impfinalfilt3",
                                "impfinalfilt4",
                                "impfinalfilt5",
                            ]
                        ):
                            if dfname in kwargs:
                                df = kwargs[dfname]
                                if df is not None:
                                    df = verifyindex(df)
                                    if not optionschanged:
                                        (
                                            reactionid_ir,
                                            reactionid_irtext,
                                            instance_ir,
                                            customrxnsmiles_ir,
                                        ) = updatereactionwidget(
                                            reactionid_ir,
                                            instance_ir,
                                            df,
                                            reactionid_irtext,
                                            customrxnsmiles_ir,
                                        )
                                        optionschanged = True
                            else:
                                kwargs[dfname] = None

            rxnview.on_click(functools.partial(on_view_clicked, view="Reaction"))
            impview.on_click(functools.partial(on_view_clicked, view="Impurity"))
            display(irmaster)

            # def traceirrxns(reactionid, instance, customrxnsmiles):
            #     global displaywidget, success

            #     reactionidui = widgets.VBox([reactionid_ir, reactionid_irtext])
            #     display(widgets.HBox([reactionidui, instance_ir, customrxnsmiles_ir]))
            #     displaywidget = True
            #     if kwargs["impfinalfilt"] is not None:
            #         kwargs["impfinalfilt"]["Relevance_morgan"] = kwargs["impfinalfilt"][
            #             "Relevance_morgan"
            #         ].apply(lambda x: round(x, 2))
            #     display(
            #         Markdown(
            #             "<h4><center><bold>12 & 13. Relevance and Conditions</bold></center></h4>"
            #         )
            #     )

            #     for i, dfname in enumerate(
            #         [
            #             "impfinalfilt",
            #             "impfinal",
            #             "analoguerxnsimpfilt",
            #             "analoguerxnsimp",
            #             "analoguerxnstempl",
            #         ]
            #     ):
            #         success = False
            #         if kwargs[dfname] is not None:
            #             df = kwargs[dfname]
            #             if (reactionid, instance) in df.index:
            #                 success = True
            #                 reactiondf = df.xs((reactionid, instance))
            #                 if i <= 4:
            #                     if "template" in reactiondf:
            #                         if isinstance(reactiondf, pd.DataFrame):
            #                             reactiondf_ = reactiondf.iloc[0]
            #                         else:
            #                             reactiondf_ = reactiondf
            #                         display(
            #                             Markdown(
            #                                 f"<h4><center><strong>9. Template</strong></center></h4>"
            #                             )
            #                         )
            #                         display(
            #                             drawReaction(
            #                                 rdChemReactions.ReactionFromSmarts(
            #                                     reactiondf_.template
            #                                 )
            #                             )
            #                         )
            #                         display(
            #                             Markdown(
            #                                 f"`Template SMILES: {reactiondf_.template}`"
            #                             )
            #                         )
            #                         display(
            #                             Markdown(
            #                                 f"<strong>Message:{reactiondf_.msg4}</strong>"
            #                             )
            #                         )
            #                 if i <= 3:
            #                     if kwargs["analoguerxnsimp"] is not None and i != 3:
            #                         reactiondf_ = kwargs["analoguerxnsimp"].xs(
            #                             (reactionid, instance)
            #                         )
            #                     else:
            #                         reactiondf_ = reactiondf
            #                     if "impurityrxn" in reactiondf_:
            #                         display(
            #                             Markdown(
            #                                 f"<h4><center><strong>10. Template Application</strong></center></h4>"
            #                             )
            #                         )
            #                         displaymsg = False
            #                         if isinstance(reactiondf_, pd.DataFrame):
            #                             for idx, reactiondf__ in reactiondf.iterrows():
            #                                 if (
            #                                     reactiondf__.impurityrxn
            #                                     and reactiondf__.impurityrxn != "Error"
            #                                 ):
            #                                     display(
            #                                         drawReaction(
            #                                             rdChemReactions.ReactionFromSmarts(
            #                                                 reactiondf__.impurityrxn,
            #                                                 useSmiles=True,
            #                                             )
            #                                         )
            #                                     )
            #                                     display(
            #                                         Markdown(
            #                                             f"`Impurity Reaction SMILES: {reactiondf__.impurityrxn}`"
            #                                         )
            #                                     )
            #                                 if not displaymsg:
            #                                     display(
            #                                         Markdown(
            #                                             f"<strong>Message:{reactiondf__.msg5}</strong>"
            #                                         )
            #                                     )
            #                                     displaymsg = True
            #                         else:
            #                             if (
            #                                 reactiondf_.impurityrxn
            #                                 and reactiondf_.impurityrxn != "Error"
            #                             ):
            #                                 for impurityrxn in reactiondf_.impurityrxn:
            #                                     display(
            #                                         drawReaction(
            #                                             rdChemReactions.ReactionFromSmarts(
            #                                                 list(impurityrxn)[0],
            #                                                 useSmiles=True,
            #                                             )
            #                                         )
            #                                     )
            #                                     display(
            #                                         Markdown(
            #                                             f"`Impurity Reaction SMILES: {list(impurityrxn)[0]}`"
            #                                         )
            #                                     )
            #                             display(
            #                                 Markdown(
            #                                     f"<strong>Message:{reactiondf_.msg5}</strong>"
            #                                 )
            #                             )
            #                     if kwargs["impfinal"] is not None:
            #                         display(
            #                             Markdown(
            #                                 f"<center><strong>{len(kwargs['impfinal'].index.unique())} reactions remaining after filtering</strong></center>"
            #                             )
            #                         )
            #                 if i <= 1:
            #                     if "msg6" in reactiondf:
            #                         display(
            #                             Markdown(
            #                                 f"<h4><center><strong>11. Impurity Cleaning</strong></center></h4>"
            #                             )
            #                         )
            #                         if isinstance(reactiondf, pd.DataFrame):
            #                             reactiondf = reactiondf.iloc[0]
            #                         if isinstance(reactiondf.msg6, pd.DataFrame):
            #                             msg6 = reactiondf.msg6.iloc[0]
            #                         else:
            #                             msg6 = reactiondf.msg6
            #                         display(
            #                             Markdown(f"<strong>Message:{msg6}</strong>")
            #                         )
            #                     if kwargs["impfinalfilt"] is not None:
            #                         display(
            #                             Markdown(
            #                                 f"<center><strong>{len(kwargs['impfinalfilt'].index.unique())} reactions remaining after filtering</strong></center>"
            #                             )
            #                         )

            #                 break

            # outip = widgets.interactive_output(
            #     traceiprxns,
            #     {
            #         "reactionid": reactionid_ip,
            #         "instance": instance_ip,
            #         "customrxnsmiles": customrxnsmiles_ip,
            #         # "analgrawcheck": analgrawcheck,
            #         # "analgcheck": analgcheck,
            #         # "analgbalcheck": analgbalcheck,
            #         # "analgmapcheck": analgmapcheck,
            #         # "analgcentcheck": analgcentcheck,
            #     },
            # )

            # display(outip)


def updatereactionwidget(
    reactionid_widget: widgets.SelectionSlider,
    instance_widget: widgets.SelectionSlider,
    df: pd.DataFrame,
    reactionidtext: widgets.Text,
    customrxnsmiles: widgets.Text,
):
    reactionids = list(df.index.get_level_values(0).unique())
    reactionid_widget.options = [""] + reactionids

    def on_update_reactionid_widget(*args):
        if reactionid_widget.value:
            instance_widget.options = df.xs(reactionid_widget.value).index.unique()
            reactionidtext.value = str(reactionid_widget.value)
            customrxnsmiles.value = ""

    def on_update_reactionidtext(*args):
        reactionid_widget.value = int(reactionidtext.value)

    def on_update_customrxnsmiles_widget(*args):
        if customrxnsmiles.value:
            reactionid_widget.value = ""
            instance_widget.value = 0

    reactionid_widget.observe(on_update_reactionid_widget, "value")
    reactionidtext.on_submit(on_update_reactionidtext)
    customrxnsmiles.on_submit(on_update_customrxnsmiles_widget)
    return reactionid_widget, reactionidtext, instance_widget, customrxnsmiles


def visreaction(rxnidentifier: str, **kwargs):
    try:
        display(
            drawReaction(
                rdChemReactions.ReactionFromSmarts(rxnidentifier, useSmiles=True)
            )
        )
    except Exception as e:
        display(Markdown("Error"))
    display(Markdown(f"`Reaction SMILES: {rxnidentifier}`"))
    if "NameDict" in kwargs and kwargs["NameDict"] is not None and kwargs["NameDict"]:
        namedict = kwargs["NameDict"]
    else:
        namedict = None
    if "LHSdata" in kwargs and kwargs["LHSdata"] is not None and kwargs["LHSdata"]:
        LHSdata = {
            specid: {
                "smiles": kwargs["LHSdata"][specid]["smiles"],
                "count": kwargs["LHSdata"][specid]["count"],
            }
            for specid in kwargs["LHSdata"]
        }
        if namedict is not None:
            for specid in LHSdata:
                LHSdata[specid].update({"name": namedict[specid]})
        display(Markdown(f"Reactant Data: {LHSdata}"))
        if "MissingReactant" in kwargs and kwargs["MissingReactant"] is not None:
            display(Markdown(f"Missing Reactants: {kwargs['MissingReactant']}"))
    if "RHSdata" in kwargs and kwargs["RHSdata"] is not None and kwargs["RHSdata"]:
        RHSdata = {
            specid: {
                "smiles": kwargs["RHSdata"][specid]["smiles"],
                "count": kwargs["RHSdata"][specid]["count"],
            }
            for specid in kwargs["RHSdata"]
        }
        if namedict is not None:
            for specid in RHSdata:
                RHSdata[specid].update({"name": namedict[specid]})
        display(Markdown(f"Product Data: {RHSdata}"))
        if "MissingProduct" in kwargs and kwargs["MissingProduct"] is not None:
            display(Markdown(f"Missing Products: {kwargs['MissingProduct']}"))
    if "Rgtdata" in kwargs and kwargs["Rgtdata"] is not None and kwargs["Rgtdata"]:
        Rgtdata = {
            specid: {
                "smiles": kwargs["Rgtdata"][specid]["smiles"],
                "count": kwargs["Rgtdata"][specid]["count"],
            }
            for specid in kwargs["Rgtdata"]
        }
        if namedict is not None:
            for specid in Rgtdata:
                Rgtdata[specid].update({"name": namedict[specid]})
        display(Markdown(f"Reagent Data: {Rgtdata}"))
        if (
            "MissingReagent" in kwargs
            and kwargs["MissingReagent"] is not None
            and kwargs["MissingReagent"]
        ):
            display(Markdown(f"Missing Reagent: {kwargs['MissingReagent']}"))
    if "Solvdata" in kwargs and kwargs["Solvdata"] is not None and kwargs["Solvdata"]:
        Solvdata = {
            specid: {
                "smiles": kwargs["Solvdata"][specid]["smiles"],
                "count": kwargs["Solvdata"][specid]["count"],
            }
            for specid in kwargs["Solvdata"]
        }
        if namedict is not None:
            for specid in Solvdata:
                Solvdata[specid].update({"name": namedict[specid]})
        display(Markdown(f"Solvent Data: {Solvdata}"))
        if (
            "MissingSolvent" in kwargs
            and kwargs["MissingSolvent"] is not None
            and kwargs["MissingSolvent"]
        ):
            display(Markdown(f"Missing Solvent: {kwargs['MissingSolvent']}"))
    if (
        "CatalystID" in kwargs
        and kwargs["CatalystID"] is not None
        and kwargs["CatalystID"]
    ):
        Catdata = kwargs["CatalystID"]
        if namedict is not None:
            Catdata = {
                catid: {"name": namedict[catid]} for catid in kwargs["CatalystID"]
            }
        display(Markdown(f"Catalyst Data: {Catdata}"))
        if (
            "MissingCatalyst" in kwargs
            and kwargs["MissingCatalyst"] is not None
            and kwargs["MissingCatalyst"]
        ):
            display(Markdown(f"Missing Catalyst: {kwargs['MissingCatalyst']}"))
    if (
        "Temperature" in kwargs
        and kwargs["Temperature"] is not None
        and kwargs["Temperature"]
    ):
        display(Markdown(f"Temperature (&deg;C): {kwargs['Temperature']}"))
    if "Pressure" in kwargs and kwargs["Pressure"] is not None and kwargs["Pressure"]:
        display(Markdown(f"Pressure (bar): {kwargs['Pressure']}"))
    if (
        "ReactionTime" in kwargs
        and kwargs["ReactionTime"] is not None
        and kwargs["ReactionTime"]
    ):
        display(Markdown(f"Reaction Time (hours): {kwargs['ReactionTime']}"))
    if "NumRefs" in kwargs and kwargs["NumRefs"] is not None and kwargs["NumRefs"]:
        display(Markdown(f"Number of References: {kwargs['NumRefs']}"))
    if "NumSteps" in kwargs and kwargs["NumSteps"] is not None and kwargs["NumSteps"]:
        display(Markdown(f"Number of Steps: {kwargs['NumSteps']}"))
    if (
        "NumStages" in kwargs
        and kwargs["NumStages"] is not None
        and kwargs["NumStages"]
    ):
        display(Markdown(f"Number of Stages: {kwargs['NumStages']}"))
    if (
        "YearPublished" in kwargs
        and kwargs["YearPublished"] is not None
        and kwargs["YearPublished"]
    ):
        display(Markdown(f"Year Published: {kwargs['YearPublished']}"))

        # if stage=='Impurity Prediction':

        # if "analoguerxns" in kwargs:
        #     analoguerxns = kwargs["analoguerxns"]
        #     analoguerxns = verifyindex(analoguerxns)

        #     def traceanaloguerxns(reactionid: int, instance: int):
        #         global success
        #         if (reactionid, instance) not in analoguerxns.index:
        #             success = False
        #             print(
        #                 f"Reaction {reactionid}, instance {instance} is not analogue"
        #             )
        #     analgout = widgets.interactive_output(
        #         traceanaloguerxns,
        #         {"reactionid": reactionid_dp, "instance": instance_dp},
        #     )

        #     if not success:
        #         display(analgout)
        #         break
        # reactiondf = analoguerxnsbal.xs((reactionid, instance))
        # if "rxnsmiles0" in reactiondf:
        #     display(
        #         Markdown(
        #             "<h3><center><strong>Reaxys Reaction</strong></center></h3>"
        #         )
        #     )
        #     rxnsmiles0 = reactiondf.rxnsmiles0
        #     display(
        #         drawReaction(
        #             rdChemReactions.ReactionFromSmarts(
        #                 rxnsmiles0, useSmiles=True
        #             )
        #         )
        #     )
        #     print(f"Reaction SMILES: {rxnsmiles0}")
        #     # if analoguerxns is not None:
        #     #     print(f"Reagents: {reactiondf.})

        # if "balrxnsmiles" in reactiondf:
        #     display(
        #         Markdown(
        #             "<h3><center><strong>Balanced Reaction</strong></center></h3>"
        #         )
        #     )
        #     balrxnsmiles = reactiondf.balrxnsmiles
        #     if balrxnsmiles != "Error":
        #         display(
        #             drawReaction(
        #                 rdChemReactions.ReactionFromSmarts(
        #                     balrxnsmiles, useSmiles=True
        #                 )
        #             )
        #         )
        #         print(f"Balanced Reaction SMILES: {balrxnsmiles}")
        #         print(f"Message: {reactiondf.msg}")
        #     else:
        #         print("Error")
        #         success = False

        #     balout = widgets.interactive_output(
        #         tracebalancerxns,
        #         {"reactionid": reactionid_dp, "instance": instance_dp},
        #     )
        #     display(balout)
        #     if not success:
        #         break
        # for i,dfname in enumerate(["analoguerxnscent","analoguerxnsassigned","analoguerxnsparsed","analoguerxnsmapped"]):

        #         def tracemappedrxns(reactionid:int,instance:int):
        #             global success
        #             display(
        #             Markdown(
        #                 "<h3><center><strong>Mapped Reaction</strong></center></h3>"
        #             )
        #             )
        #             if (reactionid,instance) not in df.index:
        #                 if i==0:

        #                     continue
        #                     print(
        #                     f"Reaction {reactionid}, instance {instance} has been removed due to an error in balancing"
        #                     )
        #                 elif i==1:

        #                 success=False
        #                 break
        #             else:
        #                 if i==0:

        #     mappedout = widgets.interactive_output(
        #         tracemappedrxns,
        #         {"reactionid": reactionid_dp, "instance": instance_dp},
        #     )
        #     display(mappedout)
        #     if not success:
        #         break
        # if 'analogue'

    # df = None
    # for dfname in [
    #     "analoguerxns",
    #     "analoguerxnsbal",
    #     "analoguerxnsmapped",
    #     "analoguerxnsparsed",
    #     "analoguerxnsassigned",
    #     "analoguerxnscent",
    #     "analoguerxnsvalid",
    #     "analoguerxnsfinal",
    #     "analoguerxnstempl",
    #     "analoguerxnsimp",
    #     "impfinal",
    #     "impfinalfilt",
    #     "impfinalfilt2",
    #     "impfinalfilt3",
    #     "impfinalfilt4",
    #     "impfinalfilt5",
    #     "summary",
    #     "summary2",
    #     "summary3",
    # ]:

    #     if dfname not in kwargs or kwargs[dfname] is None:
    #         continue
    #     else:
    #         df = verifyindex(kwargs[dfname])
    #         reactionids = list(df.index.get_level_values(0).unique())
    #         reactionid = widgets.SelectionSlider(
    #             options=reactionids, description="Reaction ID", continuous_update=False
    #         )
    #         instance = widgets.SelectionSlider(
    #             description="Instance", options=[0], continuous_update=False
    #         )

    #         def on_update_reactionid_widget(*args):
    #             instance.options = df.xs(reactionid.value).index.unique()

    #         reactionid.observe(on_update_reactionid_widget, "value")
    #         break


# def tracerxn(reactionid: int = 0, instance: int = 0,**kwargs):
#         if "analoguerxns" in kwargs:
#             analoguerxns = kwargs["analoguerxns"]
#             # analoguerxns = verifyindex(analoguerxns)
#             if (reactionid, instance) not in analoguerxns.index:
#                 print(f"Reaction {reactionid}, instance {instance} is not analogue")
#                 break
#         if "analoguerxnsbal" in kwargs:
#             analoguerxnsbal = kwargs["analoguerxnsbal"]
#             # analoguerxnsbal = verifyindex(analoguerxnsbal)
#             if (reactionid, instance) not in analoguerxnsbal.index:
#                 print(
#                     f"Reaction {reactionid}, instance {instance} does not have valid reactants or products and could not be updated"
#                 )
#                 break
#             else:
#                 reactiondf = analoguerxnsbal.xs((reactionid, instance))
#                 if "rxnsmiles0" in reactiondf:
#                     display(
#                         Markdown(
#                             "<h3><center><strong>Reaxys Reaction</strong></center></h3>"
#                         )
#                     )
#                     rxnsmiles0 = reactiondf.rxnsmiles0
#                     display(
#                         drawReaction(
#                             rdChemReactions.ReactionFromSmarts(
#                                 rxnsmiles0, useSmiles=True
#                             )
#                         )
#                     )
#                     print(f"Reaction SMILES: {rxnsmiles0}")

#                 if "balrxnsmiles" in reactiondf:
#                     display(
#                         Markdown(
#                             "<h3><center><strong>Balanced Reaction</strong></center></h3>"
#                         )
#                     )
#                     balrxnsmiles = reactiondf.balrxnsmiles
#                     if balrxnsmiles != "Error":
#                         display(
#                             drawReaction(
#                                 rdChemReactions.ReactionFromSmarts(
#                                     balrxnsmiles, useSmiles=True
#                                 )
#                             )
#                         )
#                         print(f"Balanced Reaction SMILES: {balrxnsmiles}")
#                         print(f"Message: {reactiondf.msg}")
#                     else:
#                         print("Error")

#         if "analoguerxnsmapped" in kwargs:
#             analoguerxnsmapped = kwargs["analoguerxnsmapped"]
#             # analoguerxnsmapped = verifyindex(analoguerxnsmapped)
#             display(
#                 Markdown("<h3><center><strong>Mapped Reaction</strong></center></h3>")
#             )
#             if (reactionid, instance) not in analoguerxnsmapped.index:
#                 print(
#                     f"Reaction {reactionid}, instance {instance} has been removed due to an error in balancing"
#                 )
#                 break
#             else:
#                 reactiondf = analoguerxnsmapped.xs((reactionid, instance))
#                 mappedsmiles = reactiondf.mapped_rxn
#                 if mappedsmiles != "Error":
#                     display(
#                         drawReaction(
#                             rdChemReactions.ReactionFromSmarts(
#                                 mappedsmiles, useSmiles=True
#                             )
#                         )
#                     )
#                     print(f"Mapped Reaction SMILES: {mappedsmiles}")
#                     print(f"Confidence: {reactiondf.confidence}")
#                 else:
#                     print("Error")


# if df is not None:
#     uiselect = widgets.VBox(
#         [reactionid, instance],
#         layout=widgets.Layout(border="5px solid green", width="40%"),
#     )
#     # uibutton = widgets.HBox(
#     #     [allbutton, dataminingbutton, dataprocessingbutton, impuritypredictionbutton],
#     #     layout=widgets.Layout(border="5px solid green", width="60"),
#     #     width="40%",
#     # )
#     # ui = widgets.HBox([uiselect, uibutton])
#     outmain = widgets.interactive_output(
#         tracerxn, {"reactionid": reactionid, "instance": instance}
#     )
#     # display(ui, out)
#     return uiselect, outmain


def visfragment(
    displayres=False, **kwargs
):  # Pass in either SMILES string, or a dataframe
    expandwidget = widgets.IntSlider(
        description="Expand", min=0, max=10, step=1, value=1
    )
    resformatwidget = widgets.Dropdown(
        description="Result Format",
        style={"description_width": "initial"},
        options=["smiles", "smarts"],
        value="smiles",
    )
    substanceidwidget = widgets.SelectionSlider(
        description="Species ID", options=[""], value="", continuous_update=False
    )
    substanceidtext = widgets.Text(description="Species ID", value="")
    smileswidget = widgets.SelectionSlider(
        description="SMILES", options=[""], value="", continuous_update=False
    )
    customsmileswidget = widgets.Text(
        description="Custom species SMILES",
        style={"description_width": "initial"},
    )
    fragwidget = widgets.Dropdown(description="Fragment", options=[""], value="")
    workflow = ""
    if "analoguepool" in kwargs and isinstance(kwargs["analoguepool"], pd.DataFrame):
        analoguepool = kwargs["analoguepool"]
        analoguepool = verifyindex(analoguepool, idxcol=["SubstanceID"])
        substanceidwidget.options = [""] + list(analoguepool.index)
        smileswidget.options = [""] + list(analoguepool.Smiles.values.unique())
        workflow = "analoguepool"
    else:
        analoguepool = None
    if "fragdbsource" in kwargs and isinstance(kwargs["fragdbsource"], pd.DataFrame):
        fragdb = kwargs["fragdbsource"]
        fragdb = verifyindex(fragdb, idxcol=["FragmentSmiles", "SubstanceID"])
        if not workflow:
            substanceidwidget.options = [""] + list(
                fragdb.index.get_level_values(1).unique()
            )
            smileswidget.options = [""] + list(fragdb.Smiles.values.unique())
            workflow = "fragdb"
    else:
        fragdb = None
    if "substancesource" in kwargs and isinstance(
        kwargs["substancesource"], pd.DataFrame
    ):
        substancedb = kwargs["substancesource"]
        substancedb = verifyindex(substancedb, idxcol="SubstanceID")
        if not workflow:
            substanceidwidget.options = [""] + list(substancedb.index)
            smileswidget.options = [""] + list(substancedb.Smiles.values.unique())
            workflow = "substancedb"
    else:
        substancedb = None
    if "analoguefrag" in kwargs:
        analoguefrag = kwargs["analoguefrag"]
    else:
        analoguefrag = None

    def on_update_substance_id(*args):
        if not substanceidwidget.value:
            substanceidtext.value = ""
            smileswidget.value = ""
            fragwidget.value = ""
            fragwidget.options = [""]

        else:
            substanceidtext.value = str(substanceidwidget.value)
            if workflow == "analoguepool":
                smiles = analoguepool.loc[substanceidwidget.value].Smiles
                if isinstance(smiles, pd.Series):
                    smiles = smiles.values[0]
                smileswidget.value = smiles
            elif workflow == "fragdb":
                smiles = fragdb.xs(substanceidwidget.value, level=1).Smiles
                if isinstance(smiles, pd.Series):
                    smiles = smiles.values[0]
                smileswidget.value = smiles
            else:
                smiles = substancedb.loc[substanceidwidget.value].Smiles
                if isinstance(smiles, pd.Series):
                    smiles = smiles.values[0]
                smileswidget.value = smiles
            if (
                fragdb is not None
                and expandwidget.value == 1
                and resformatwidget.value == "smiles"
            ):
                fraglist = list(fragdb.xs(substanceidwidget.value, level=1).index)
            else:  # Manual
                if "." in smileswidget.value:  # Mixture
                    frags = getmixturefrags(
                        smileswidget.value,
                        expand=expandwidget.value,
                        resFormat=resformatwidget.value,
                    )
                else:
                    frags = getCarrierFrags0(
                        smileswidget.value,
                        expand=expandwidget.value,
                        resFormat=resformatwidget.value,
                    )
                fraglist = list(Counter(frags).keys())
            fragwidget.options = [""] + fraglist
            if analoguefrag is not None and analoguefrag in fragwidget.options:
                fragwidget.value = analoguefrag

    def on_update_smiles(*args):
        if not smileswidget.value:
            substanceidwidget.value = ""
            fragwidget.value = ""
            fragwidget.options = [""]
        else:
            if workflow == "analoguepool":
                substanceidwidget.value = analoguepool.loc[
                    analoguepool.Smiles == smileswidget.value
                ].index[0]
            elif workflow == "fragdb":
                substanceidwidget.value = fragdb.loc[
                    fragdb.Smiles == smileswidget.value
                ].index.get_level_values(1)[0]
            else:
                substanceidwidget.value = substancedb.loc[
                    substancedb.Smiles == smileswidget.value
                ].index[0]
            if (
                fragdb is not None
                and expandwidget.value == 1
                and resformatwidget.value == "smiles"
            ):
                fraglist = list(fragdb.xs(substanceidwidget.value, level=1).index)
            else:  # Manual
                if "." in smileswidget.value:  # Mixture
                    frags = getmixturefrags(
                        smileswidget.value,
                        expand=expandwidget.value,
                        resFormat=resformatwidget.value,
                    )
                else:
                    frags = getCarrierFrags0(
                        smileswidget.value,
                        expand=expandwidget.value,
                        resFormat=resformatwidget.value,
                    )
                fraglist = list(Counter(frags).keys())
            fragwidget.options = [""] + fraglist
            if analoguefrag is not None and analoguefrag in fragwidget.options:
                fragwidget.value = analoguefrag

    def on_update_substanceidtext(*args):
        substanceidwidget.value = int(substanceidtext.value)

    def on_update_custom_smiles(*args):
        present = False
        try:
            customsmiles = Chem.MolToSmiles(molfromsmiles(customsmileswidget.value))
            # raise Exception
            if workflow == "analoguepool":
                substanceidwidget.value = analoguepool.loc[
                    analoguepool.Smiles == customsmiles
                ].index[0]
            elif workflow == "fragdb":
                substanceidwidget.value = fragdb.loc[
                    fragdb.Smiles == customsmiles
                ].index.get_level_values(1)[0]
            else:
                substanceidwidget.value = substancedb.loc[
                    substancedb.Smiles == customsmiles
                ].index[0]
            customsmileswidget.value = customsmiles
            smileswidget.value = customsmileswidget.value
            present = True
            if fragdb is not None:
                fraglist = list(fragdb.xs(substanceidwidget.value, level=1).index)
            else:
                raise Exception
        except Exception as e:  # Custom smiles does not exist
            if not present:
                substanceidwidget.value = ""
                smileswidget.value = ""
            if "." in customsmileswidget.value:  # Mixture
                frags = getmixturefrags(
                    customsmileswidget.value,
                    expand=expandwidget.value,
                    resFormat=resformatwidget.value,
                )
            else:
                frags = getCarrierFrags0(
                    customsmileswidget.value,
                    expand=expandwidget.value,
                    resFormat=resformatwidget.value,
                )
            fraglist = list(Counter(frags).keys())
            fragwidget.options = [""] + fraglist
            if analoguefrag is not None and analoguefrag in fragwidget.options:
                fragwidget.value = analoguefrag

    def on_update_expand(*args):
        smiles = smileswidget.value
        if not smiles:
            if customsmileswidget.value:
                smiles = customsmileswidget.value
        if smiles:
            if (
                expandwidget.value != 1
                or resformatwidget.value != "smiles"
                or fragdb is None
                or (not substanceidwidget.value)
            ):
                if "." in smiles:  # Mixture
                    frags = getmixturefrags(
                        smiles,
                        expand=expandwidget.value,
                        resFormat=resformatwidget.value,
                    )
                else:
                    frags = getCarrierFrags0(
                        smiles,
                        expand=expandwidget.value,
                        resFormat=resformatwidget.value,
                    )
                fraglist = list(Counter(frags).keys())
            else:
                fraglist = list(fragdb.xs(substanceidwidget.value, level=1).index)
            fragwidget.options = [""] + fraglist
            if analoguefrag is not None and analoguefrag in fragwidget.options:
                fragwidget.value = analoguefrag

    def on_update_resformat(*args):
        smiles = smileswidget.value
        if not smiles:
            if customsmileswidget.value:
                smiles = customsmileswidget.value
        if smiles:
            if (
                resformatwidget.value != "smiles"
                or expandwidget.value != 1
                or fragdb is None
                or (not substanceidwidget.value)
            ):
                if "." in smiles:  # Mixture
                    frags = getmixturefrags(
                        smiles,
                        expand=expandwidget.value,
                        resFormat=resformatwidget.value,
                    )
                else:
                    frags = getCarrierFrags0(
                        smiles,
                        expand=expandwidget.value,
                        resFormat=resformatwidget.value,
                    )
                fraglist = list(Counter(frags).keys())
            else:
                fraglist = list(fragdb.xs(substanceidwidget.value, level=1).index)
            fragwidget.options = [""] + fraglist
            if analoguefrag is not None and analoguefrag in fragwidget.options:
                fragwidget.value = analoguefrag

    # smileslabel=Label('Input a substance SMILES here:')
    # smileswidget=widgets.HBox([smileslabel,widgets.Text(description="SMILES")])
    def f(frag: str, spec: str, customspec: str):
        if frag:
            if spec:
                img, count = highlightsubstruct(spec, [frag], returncount=True)
            elif customspec:
                img, count = highlightsubstruct(customspec, [frag], returncount=True)
            display(img)
            if spec:
                print(f"Species SMILES: {spec}")
            elif customspec:
                print(f"Species SMILES: {customspec}")
            print(f"Fragment identifier: {frag}")
            print(f"Fragment count: {count[0]}")

    substanceidwidget.observe(on_update_substance_id, "value")
    substanceidtext.on_submit(on_update_substanceidtext)
    smileswidget.observe(on_update_smiles, "value")
    customsmileswidget.observe(on_update_custom_smiles, "value")
    expandwidget.observe(on_update_expand, "value")
    resformatwidget.observe(on_update_resformat, "value")
    out = widgets.interactive_output(
        f,
        {
            "frag": fragwidget,
            "spec": smileswidget,
            "customspec": customsmileswidget,
        },
    )
    substanceidui = widgets.VBox(
        [substanceidwidget, substanceidtext],
        layout=widgets.Layout(width="50%", height="auto"),
    )
    smilesui = widgets.VBox(
        [smileswidget, customsmileswidget],
        layout=widgets.Layout(width="50%", height="auto"),
    )
    combinedui = widgets.HBox(
        [substanceidui, smilesui],
        layout=widgets.Layout(border="5px solid green", width="100%", height="auto"),
    )
    fragui = widgets.HBox(
        [fragwidget, expandwidget, resformatwidget],
        layout=widgets.Layout(border="5px solid green", width="100%", height="auto"),
    )
    # masterfragui = widgets.VBox([combinedui, fragui])
    if displayres:
        display(Markdown("<h3><center><strong>Species</strong></center></h3>"))
        display(combinedui)
        display(Markdown("<h3><center><strong>Fragments</strong></center></h3>"))
        display(fragui)
        display(out)
    else:
        return combinedui, fragui, out


# def visfragment(
#     userinput: Union[pd.DataFrame, str, Dict],
#     expand: int = 1,
#     resFormat: str = "smiles",
#     **kwargs,
# ):

#     fragdb = False
#     substancedb = False
#     if isinstance(userinput, pd.DataFrame):  # User has passed in a dataframe
#         for i, identifier in enumerate(["FragmentSmiles", "FragmentSmarts", "Smiles"]):
#             if (
#                 identifier in userinput.columns
#                 or identifier in userinput.index.names
#                 or userinput.index.name == identifier
#             ):
#                 if i < 2:
#                     fragdb = True
#                 elif ">1 Compound" in userinput.columns:
#                     substancedb = True
#                 if (
#                     identifier in userinput.index.names
#                     or userinput.index.name == identifier
#                 ):
#                     userinput.reset_index(inplace=True)
#                 break
#         if not fragdb and not substancedb:
#             raise CustomError("User input dataframe could not be processed")
#         if "SubstanceID" not in userinput.columns:
#             userinput.reset_index(inplace=True)
#         substanceids = userinput.SubstanceID.unique()
#     else:
#         identifier = ""
#         substanceids = [0]

#     a = widgets.SelectionSlider(description="SubstanceID", options=substanceids)

#     def f(a):
#         substanceid = a
#         if fragdb:
#             fraginfo = userinput.loc[userinput.SubstanceID == substanceid][
#                 [identifier, "count"]
#             ].to_dict(orient="split")["data"]
#             smiles = userinput.loc[userinput.SubstanceID == substanceid].Smiles.iloc[0]
#         else:
#             if substancedb:
#                 smiles = userinput.loc[userinput.SubstanceID == substanceid][
#                     identifier
#                 ].iloc[0]
#             else:
#                 smiles = userinput
#             try:
#                 if "." in smiles:  # Mixture
#                     fraglist = getmixturefrags(
#                         smiles, expand=expand, resFormat=resFormat
#                     )
#                 else:
#                     fraglist = getCarrierFrags0(
#                         smiles, expand=expand, resFormat=resFormat
#                     )
#             except Exception as e:
#                 raise CustomError("User input could not be processed")
#             fraginfo = Counter(fraglist).most_common()
#         for frag, count in fraginfo:
#             print(f"Species ID: {substanceid}")
#             print(f"Species SMILES: {smiles}")
#             print(f"Fragment: {frag}")
#             print(f"Count: {count}")
#             display(highlightsubstruct(smiles, pattlist=[frag]))

#     out = widgets.interactive_output(f, {"a": a})
#     display(a)
#     display(out)


def visoutput(analoguerxns: pd.DataFrame):
    if analoguerxns.index.name or any(analoguerxns.index.names):
        analoguerxns.reset_index(inplace=True)
    a = widgets.IntSlider(min=0, max=len(analoguerxns) - 1, description="Row Number")

    def f(a):
        reaxysID = analoguerxns.iloc[a].ReactionID
        print("Reaxys reaction " + str(reaxysID) + ":")
        if "rxnsmiles0" in analoguerxns.dtypes:
            rxnsmiles0 = analoguerxns.iloc[a].rxnsmiles0
            display(
                drawReaction(
                    rdChemReactions.ReactionFromSmarts(rxnsmiles0, useSmiles=True)
                )
            )
        if "balrxnsmiles" in analoguerxns.dtypes:
            balrxnsmiles = analoguerxns.iloc[a].balrxnsmiles
            msg = analoguerxns.iloc[a].msg
            print("Balancing algorithm output: " + msg)
            print("Balanced reaction:")
            if balrxnsmiles != "Error":
                display(
                    drawReaction(
                        rdChemReactions.ReactionFromSmarts(balrxnsmiles, useSmiles=True)
                    )
                )
                print(balrxnsmiles)
            else:
                print("Error")
        if "mapped_rxn" in analoguerxns.dtypes:
            mappedrxn = analoguerxns.iloc[a].mapped_rxn
        else:
            try:
                mappedrxn = maprxn([balrxnsmiles])[0]["mapped_rxn"]
            except Exception:
                mappedrxn = "Error"
        print("Mapped reaction:")
        if mappedrxn != "Error":
            display(
                drawReaction(
                    rdChemReactions.ReactionFromSmarts(mappedrxn, useSmiles=True)
                )
            )
        else:
            print(mappedrxn)
        if "msg1" in analoguerxns.dtypes:
            msg1 = analoguerxns.iloc[a].msg1
            print("Mapping validity: " + msg1)
        if "template" in analoguerxns.dtypes:
            template = analoguerxns.iloc[a].template
            msg4 = analoguerxns.iloc[a].msg4
            print("Template message: " + msg4)
            print("Template:")
            display(
                drawReaction(
                    rdChemReactions.ReactionFromSmarts(template, useSmiles=True)
                )
            )
            print("Template SMARTS: " + template)

    out = widgets.interactive_output(f, {"a": a})
    display(a)
    display(out)


def visoutput2(analoguerxns: pd.DataFrame):
    #     breakpoint()
    if analoguerxns.index.name or any(analoguerxns.index.names):
        analoguerxns.reset_index(inplace=True)
    b = widgets.IntSlider(
        min=min(analoguerxns.ReactionID),
        max=max(analoguerxns.ReactionID),
        description="Reaxys ID",
    )
    c = widgets.IntSlider(
        min=0,
        max=len(analoguerxns.loc[analoguerxns.ReactionID == b]),
        description="Instance",
    )

    def f(b, c):
        reaxysID = b
        inst = c
        print("Reaxys reaction " + str(reaxysID) + ":")
        if "rxnsmiles0" in analoguerxns.dtypes:
            rxnsmiles0 = (
                analoguerxns.loc[analoguerxns.ReactionID == b].iloc[c].rxnsmiles0
            )
            display(
                drawReaction(
                    rdChemReactions.ReactionFromSmarts(rxnsmiles0, useSmiles=True)
                )
            )
        if "balrxnsmiles" in analoguerxns.dtypes:
            balrxnsmiles = (
                analoguerxns.loc[analoguerxns.ReactionID == b].iloc[c].balrxnsmiles
            )
            msg = analoguerxns.loc[analoguerxns.ReactionID == b].iloc[c].msg
            print("Balancing algorithm output: " + msg)
            print("Balanced reaction:")
            if balrxnsmiles != "Error":
                display(
                    drawReaction(
                        rdChemReactions.ReactionFromSmarts(balrxnsmiles, useSmiles=True)
                    )
                )
                print(balrxnsmiles)
            else:
                print("Error")
        if "mapped_rxn" in analoguerxns.dtypes:
            mappedrxn = (
                analoguerxns.loc[analoguerxns.ReactionID == b].iloc[c].mapped_rxn
            )
        else:
            try:
                mappedrxn = maprxn([balrxnsmiles])[0]["mapped_rxn"]
            except Exception:
                mappedrxn = "Error"
        print("Mapped reaction:")
        if mappedrxn != "Error":
            display(
                drawReaction(
                    rdChemReactions.ReactionFromSmarts(mappedrxn, useSmiles=True)
                )
            )
        else:
            print(mappedrxn)
        if "msg1" in analoguerxns.dtypes:
            msg1 = analoguerxns.loc[analoguerxns.ReactionID == b].iloc[c].msg1
            print("Mapping validity: " + msg1)
        if "template" in analoguerxns.dtypes:
            template = analoguerxns.loc[analoguerxns.ReactionID == b].iloc[c].template
            msg4 = analoguerxns.loc[analoguerxns.ReactionID == b].iloc[c].msg4
            print("Template message: " + msg4)
            print("Template:")
            display(
                drawReaction(
                    rdChemReactions.ReactionFromSmarts(template, useSmiles=True)
                )
            )

    out = widgets.interactive_output(f, {"b": b, "c": c})
    display(b)
    display(c)
    display(out)


# demo = pd.read_pickle(
#     "/home/aa2133/Impurity-Project/Reaxys_Data/SubstanceSmiles.pickle"
# )
# visfragment(demo[:10])

# def balance_rxn(rxn: str, **kwargs):
#     """_summary_

#     Args:
#         rxn (str): _description_
#         **kwargs: List of possible variables indicated within function scope.
#     """
#     IP = {
#         "refbalrxns": None,  # If prior results need to be included
#         "reagents": [],  # Only if one reaction is inputted
#         "solvents": [],  # Only if one reaction is inputted
#         "coefflim": 6,  # Maximum tolerable stoichiometric coefficient
#         "reaxys_update": True,  # If Reaxys has been used to update
#         "includesolv": True,  # If solvents need to be used
#         "usemapper": True,  # If mapper needs to be used
#         "helpprod": True,  # If help products need to be used
#         "helpreact": False,  # If help reactants need to be used
#         "addrctonly": False,  # If only reactants should be included for balancing
#         "ignoreH": False,  # If hydrogens are to be ignored when balancing
#         "ncpus": 1,  # Number of CPUs for computation
#         "restart": True,  # Restart distributed cluster,
#         "helpprod": True,  # Use help products or not
#         "helpreact": False,  # Use help reactants or not
#         "hc_prod": hc_Dict,  # Help compound dictionary,
#         "hc_react": None,
#         "first": False,
#     }  # Help reactant dictionary
#     IP = {**IP, **kwargs}
#     Rdata, Pdata, Rgtdata, Solvdata = getspecdat_rxn(
#         rxn, reagents=IP["reagents"], solvents=IP["solvents"]
#     )
#     rxn = pd.DataFrame(
#         [
#             {
#                 "ReactionID": 0,
#                 "Instance": 0,
#                 "NumSteps": 1,
#                 "NumStages": 1,
#                 "NumRefs": 1,
#                 "Rdata": Rdata,
#                 "Pdata": Pdata,
#                 "Rgtdata": Rgtdata,
#                 "Solvdata": Solvdata,
#                 "hc_prod": IP["hc_prod"],
#                 "hc_react": IP["hc_react"],
#             }
#         ]
#     )
#     IP["addedspecies"] = [i for i in Rdata]

#     balrxnsraw, balancedrxns = balance_analogue_(rxn, **IP)

#     mappedrxns=map_rxns(balancedrxns,ncpus=IP['ncpus'],restart=IP['restart'],reaxys_update=IP['reaxys_update'])
#     checkedrxns=checkrxns(mappedrxns,reaxys_update=IP['reaxys_update'],ncpus=IP['ncpus'])
#     changedrxns=checkedrxns.loc[(checkedrxns.msg1.str.contains('Unmapped')) | (checkedrxns.msg1.str.contains('unmapped'))]
#     changedrxns=changedrxns.loc[~(changedrxns.msg1.str.contains('Error')) & ~(changedrxns.msg1.str.contains('Invalid')) & ~(changedrxns.msg1.str.contains('Mandatory',case=False,na=False)) & ~(changedrxns.msg1.str.contains('discrepancy',case=False,na=False))]
#     if not changedrxns.empty:
#         changedrxns=updaterxns(changedrxns,hc_prod=IP['hc_prod'],analoguerxns=rxn,ncpus=IP['ncpus'])
#         checkedrxns.update(changedrxns[['mapped_rxn','confidence','balrxnsmiles','msg','LHS','RHS','hcrct','hcprod','LHSdata','RHSdata','msg1']])
#     return checkedrxns

# For a custom reaction SMILES string
