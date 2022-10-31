# %load ReaxysAPIv2.py

## Using Reaxys script

# -*- coding: utf-8 -*-
# Python wrapper for the Reaxys API
#
# Version: 1.1.0-beta.2
#
# Author:  Dr. Sebastian Radestock, Elsevier
# Author:  Dr. Alexander Riemer, Elsevier
# Author:  Dr. Markus Fischer, Elsevier
# Date:    July 26th, 2019
# Change Log 1.1.0-beta.1, July 26th, 2019
# A. Support for Python 3
# B. get_field_content modifications
# B.1. returns values for elements with highlights
# B.2. new method argument highlight_only. If True will return only a value if field contains highlights
#
# Change Log 1.1.0-beta.2, July 26th, 2019
# A. Method retrieve now supports clustering hitsets
# A.1. Added optional arguments dbname and context, that are required to formulate group by statements

import argparse
import http.cookiejar
import os
import os.path
import re
import string
import sys

# Newly added
import time
import traceback
import xml.dom.minidom
from urllib.request import Request, urlopen

import ray
from lxml import etree

#%% Parallel set-up (Comment if not desired) #%%
# def initray(restart=True,num_cpus=16,log_to_driver=False):
#     '''
#     Initiates cluster of CPUs

#     '''
#     if restart:
#         ray.shutdown()
#     ray.init(num_cpus=num_cpus,log_to_driver=log_to_driver)

# initray(num_cpus=5) #Creating a cluster of 5 CPUs


@ray.remote(
    num_cpus=1
)  # Each actor can only use 1 CPU. Comment line if parallel execution is not desired.
class Reaxys_API:
    def __init__(self, proxy=None, port=None):

        self.url = ""
        self.headers = {"Content-type": 'text/xml; charset="UTF-8"'}
        self.callername = ""
        self.sessionid = ""
        self.resultname = ""
        self.resultsize = ""
        self.citationset = ""
        self.citationcount = ""
        self.proxy = proxy
        self.port = port
        self.active = False

        # Set True for verbose output:
        self.debug = False

    def __del__(self):
        if self.active == True:
            print(f"INFO: Closing session {self.sessionid}")
            self.disconnect()

    def _get_resultname(self, response_xml):

        response_dom = xml.dom.minidom.parseString(response_xml)

        # Length of response_dom.getElementsByTagName("resultsname") should always be 1.
        # Node resultsname should not conatin subnodes.
        try:
            resultname = (
                response_dom.getElementsByTagName("resultname")[0]
                .childNodes[0]
                .nodeValue
            )
        except IndexError:
            resultname = None
        return resultname

    def _get_resultsize(self, response_xml):

        response_dom = xml.dom.minidom.parseString(response_xml)

        # Length of response_dom.getElementsByTagName("resultsize") should always be 1.
        # Node resultsize should not conatin subnodes.
        try:
            resultsize = (
                response_dom.getElementsByTagName("resultsize")[0]
                .childNodes[0]
                .nodeValue
            )
        except IndexError:
            resultsize = None

        return resultsize

    def _get_citationset(self, response_xml):

        response_dom = xml.dom.minidom.parseString(response_xml)

        # Length of response_dom.getElementsByTagName("citationset") should always be 1.
        # Node citationset should not conatin subnodes.
        return (
            response_dom.getElementsByTagName("citationset")[0].childNodes[0].nodeValue
        )

    def _get_citationcount(self, response_xml):

        response_dom = xml.dom.minidom.parseString(response_xml)

        # Length of response_dom.getElementsByTagName("citationcount") should always be 1.
        # Node citationcount should not conatin subnodes.
        return (
            response_dom.getElementsByTagName("citationcount")[0]
            .childNodes[0]
            .nodeValue
        )

    def get_facts_availability(self, response_xml, field):

        facts_availability = "0"

        response_dom = xml.dom.minidom.parseString(response_xml)

        facts = response_dom.getElementsByTagName("facts")[0]
        for fact in facts.childNodes:
            if 'name="' + field + '"' in fact.toxml():
                facts_availability = fact.childNodes[0].nodeValue.split("(")[0]

        return facts_availability

    def get_field_content(
        self, response_xml, field, highlight_only=False, plaintext=False
    ):

        field_content = []

        response_dom = xml.dom.minidom.parseString(response_xml)

        for element in response_dom.getElementsByTagName(field):
            if plaintext:
                field_content.append(
                    "".join([child.toxml() for child in element.childNodes])
                )

            # Concatenate text values if highlight is present
            elif element.getAttribute("highlight") == "true":
                field_content.append(
                    "".join(
                        [
                            e.data
                            if type(e) == xml.dom.minidom.Text
                            else e.childNodes[0].data
                            for e in element.childNodes
                        ]
                    )
                )

            # If node contains further sub-nodes: return full xml.
            elif len(element.childNodes) > 1 and highlight_only is False:
                field_content.append(element.toxml())

            # If node does not conatin further sub-nodes: return node value.
            elif len(element.childNodes) == 1 and highlight_only is False:
                field_content.append(element.childNodes[0].nodeValue)

        return field_content

    def connect(self, url, username, password, callername):  # url_main,

        self.url = url
        self.callername = callername
        cookies = http.cookiejar.CookieJar()

        connect_template = """<?xml version="1.0"?>
          <xf>
            <request caller="%s">
              <statement command="connect" username="%s" password="%s"/>
            </request>
          </xf>\n"""
        payload = connect_template % (callername, username, password)
        data = payload.encode()

        # Header reset.
        self.headers = {"Content-type": 'text/xml; charset="UTF-8"'}

        # ELSAPI support
        self.headers["X-ELS-APIKey"] = callername
        self.headers["Accept"] = "*/*"
        request = Request(self.url, data=data, headers=self.headers)

        if self.debug:
            print("-----------------------\nQuery headers from connect:")
            print(self.headers)
            print("-----------------------\nQuery from connect:")
            print(payload)

        response = urlopen(request)
        response_xml = response.read()

        if self.debug:
            print("-----------------------\nResponse headers from connect:")
            print(response.info())
            print("-----------------------\nResponse from connect:")
            print(response_xml)

        # Get sessionid.
        response_dom = xml.dom.minidom.parseString(response_xml)
        element = response_dom.getElementsByTagName("sessionid")
        self.sessionid = element[0].childNodes[0].nodeValue

        # Cookies are read from the response and stored in self.header
        #     which is used as a request header for subsequent requests.
        cookies.extract_cookies(response, request)

        # Cookie handling 3.0: Simply store and resend ALL cookies received from server
        self.headers["Cookie"] = "; ".join(
            re.findall(r"(?<=Cookie ).*?=\S*", str(cookies))
        )
        self.active = True

    def disconnect(self):

        disconnect_template = """<?xml version="1.0"?>
          <xf>
            <request caller="%s">
              <statement command="disconnect" sessionid="%s"/>
            </request>
          </xf>\n"""
        payload = disconnect_template % (self.callername, self.sessionid)
        data = payload.encode()

        request = Request(self.url, data=data, headers=self.headers)

        if self.debug:
            print("-----------------------\nQuery headers from disconnect:")
            print(self.headers)
            print("-----------------------\nQuery from disconnect:")
            print(payload)

        response = urlopen(request)
        response_xml = response.read()

        if self.debug:
            print("-----------------------\nResponse headers from disconnect:")
            print(response.info())
            print("-----------------------\nResponse from disconnect:")
            print(response_xml)

        self.active = False

    def select(self, dbname, context, where_clause, order_by, options):
        #         breakpoint()
        select_template = """<?xml version="1.0" encoding="UTF-8"?>
          <xf>
            <request caller="%s" sessionid="">
              <statement command="select"/>
              <select_list>
                <select_item/>
              </select_list>
              <from_clause dbname="%s" context="%s">
              </from_clause>
              <where_clause>%s</where_clause>
              <order_by_clause>%s</order_by_clause>
              <options>%s</options>
            </request>
          </xf>\n"""
        payload = select_template % (
            self.callername,
            dbname,
            context,
            where_clause,
            order_by,
            options,
        )
        data = payload.encode()
        request = Request(self.url, data=data, headers=self.headers)

        if self.debug:
            print("-----------------------\nQuery headers from select:")
            print(self.headers)
            print("-----------------------\nQuery from select:")
            print(payload)

        response = urlopen(request)
        response_xml = response.read()

        if self.debug:
            print("-----------------------\nResponse headers from select:")
            print(response.info())
            print("-----------------------\nResponse from select:")
            print(response_xml)

        self.resultname = self._get_resultname(response_xml)
        self.resultsize = self._get_resultsize(response_xml)

        if ("NO_CORESULT" not in options) and ("C" not in context):
            self.citationset = self._get_citationset(response_xml)
            self.citationcount = self._get_citationcount(response_xml)

        return response_xml

    def expand(self, dbname, first_item, last_item, where_clause):

        select_template = """<?xml version="1.0" encoding="UTF-8"?>
          <xf>
            <request caller="%s" sessionid="%s">
              <statement command="expand"/>
              <from_clause dbname="%s" first_item="%s" last_item="%s">
              </from_clause>
              <where_clause>%s</where_clause>
            </request>
          </xf>\n"""
        payload = select_template % (
            self.callername,
            self.sessionid,
            dbname,
            first_item,
            last_item,
            where_clause,
        )
        data = payload.encode()
        request = Request(self.url, data=data, headers=self.headers)

        if self.debug:
            print("-----------------------\nQuery headers from expand:")
            print(self.headers)
            print("-----------------------\nQuery from expand:")
            print(payload)

        response = urlopen(request)
        response_xml = response.read()

        if self.debug:
            print("-----------------------\nResponse headers from expand:")
            print(response.info())
            print("-----------------------\nResponse from expand:")
            print(response_xml)

        return response_xml

    def post(self, payload):

        data = payload.encode()
        request = Request(self.url, data=data, headers=self.headers)

        if self.debug:
            print("-----------------------\nQuery headers from post:")
            print(self.headers)
            print("-----------------------\nQuery from post:")
            print(payload)

        response = urlopen(request)
        response_xml = response.read()

        if self.debug:
            print("-----------------------\nResponse headers from post:")
            print(response.info())
            print("-----------------------\nResponse from post:")
            print(response_xml)

    def retrieve(
        self,
        resultname,
        select_items,
        first_item,
        last_item,
        order_by,
        group_by,
        group_item,
        options,
        dbname=None,
        context=None,
    ):
        # if group_by is given, please provide group_item, e.g. "1" or "1,2"

        if group_by != "":
            grouplist = 'grouplist="' + group_item + '"'
        else:
            grouplist = ""

        db_template = ""
        if dbname is not None:
            db_template = 'dbname="%s"' % dbname

        context_template = ""
        if context is not None:
            context_template = 'context="%s"' % context

        select_item_template = """                <select_item>%s</select_item>\n"""
        select_template = """<?xml version="1.0" encoding="UTF-8"?>
          <xf>
            <request caller="%s" sessionid="%s">
              <statement command="select"/>
              <select_list>\n"""
        for index in range(0, len(select_items)):
            select_template = select_template + select_item_template % (
                select_items[index]
            )
        select_template = (
            select_template
            + """              </select_list>
              <from_clause %s %s resultname="%s" %s first_item="%s" last_item="%s">
              </from_clause>
              <order_by_clause>%s</order_by_clause>
              <group_by_clause>%s</group_by_clause>
              <options>%s</options>
            </request>
          </xf>\n"""
        )
        payload = select_template % (
            self.callername,
            self.sessionid,
            db_template,
            context_template,
            resultname,
            grouplist,
            first_item,
            last_item,
            order_by,
            group_by,
            options,
        )
        data = payload.encode()

        request = Request(self.url, data=data, headers=self.headers)

        if self.debug:
            print("-----------------------\nQuery headers from retrieve:")
            print(self.headers)
            print("-----------------------\nQuery from retrieve:")
            print(payload)

        response = urlopen(request)
        response_xml = response.read().decode()

        if self.debug:
            print("-----------------------\nResponse headers from retrieve:")
            print(response.info())
            print("-----------------------\nResponse from retrieve:")
            print(response_xml)

        return response_xml

    #%%%%%%%%%%%%%% Added batch reaction extraction functionality #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def get_species(self, xml, tag, tagname, parenttag=None):
        """
        Retrieve species ID, any missing information (eg./ plaintext without ID), species named
        eg./ get_species(xml=rxd,tag='RXD.RGTXRN',tagname='RXD.RGT',parenttag='RXD03')

        """
        missingspec = []
        spec = []
        specnames = []
        if parenttag is not None:
            for i, specgroup in enumerate(
                list(self.get_field_content(xml, parenttag)), start=1
            ):
                specl = list(self.get_field_content(specgroup, tag))
                speclname = list(
                    self.get_field_content(specgroup, tagname, plaintext=True)
                )
                if speclname and not specl:  # Name but no ID
                    missingspec += speclname
                else:
                    #                 spec+=[int(ID) for ID in specl]
                    spec += specl
                    if speclname:
                        specnames += speclname
                    else:
                        specnames += [None]  # ID but no name
        else:
            specl = list(self.get_field_content(xml, tag))
            speclname = list(self.get_field_content(xml, tagname, plaintext=True))
            if speclname and not specl:  # Name but no ID
                missingspec += speclname
            else:
                #             spec+=[int(ID) for ID in specl]
                spec += specl
                if speclname:
                    specnames += speclname
                else:
                    specnames += [None]  # ID but no name
        return spec, missingspec, specnames

    def get_conds(self, xml, tag):
        cond = list(self.get_field_content(xml, tag, plaintext=True))
        return cond

    def get_yield(self, xml):
        yieldinfo = {}
        yieldgroups = list(self.get_field_content(xml, "RXD01"))
        if not yieldgroups:
            return yieldinfo
        for yieldgroup in yieldgroups:
            yieldcompd = list(self.get_field_content(yieldgroup, "RXD.YXRN"))
            yieldval = list(self.get_field_content(yieldgroup, "RXD.NYD"))
            yieldinfo.update(dict(zip(yieldcompd, yieldval)))
        return yieldinfo

    def get_details(self, reaction, coreinfo, current_rindex, rxd_chunksize):
        """
        Process variations for a given reaction and collect facts as defined in the function get_facts.
        If more than 20 (rxd_chunksize) variations exist for a reaction, the additional variations are retrieved
        and processed.

        :param reaction: string representation of the reaction XML as returned by the API
        :param current_rindex: the position of the reaction within the Reaxys hit set.
                               Required if more than 20 (rxd_chunksize) have to be retrieved
        :return: list of dicts where each dict represents the details for one variation
        coreinfo: List (ReactantID,ProductID,ReactantNames,ProductNames)
        """
        #     breakpoint()
        #     print("Current Reaction Index %s" % (current_rindex,))
        details = []
        rxd_list = list(self.get_field_content(reaction, "RXD"))
        # Citations and conditions (e.g reagents) are associated with
        # variations (RXD). Therefore we must loop over all variations to retrieve
        # all conditions and associated DOIs.
        #     print(f"Processing RXDs from 1 to {rxd_chunksize}")
        for rxd in rxd_list:
            numsteps = self.get_conds(rxd, "RXD.STP")[0]
            rgt, missingrgt, rgtnames = self.get_species(
                rxd, "RXD.RGTXRN", "RXD.RGT", parenttag="RXD03"
            )
            cat, missingcat, catnames = self.get_species(
                rxd, "RXD.CATXRN", "RXD.CAT", parenttag="RXD04"
            )
            solv, missingsolv, solvnames = self.get_species(
                rxd, "RXD.SOLXRN", "RXD.SOL", parenttag="RXD05"
            )
            temp = self.get_conds(rxd, "RXD.T")
            press = self.get_conds(rxd, "RXD.P")
            time = self.get_conds(rxd, "RXD.TIM")
            condnotes = ", ".join(
                self.get_conds(rxd, "RXD.COND") + self.get_conds(rxd, "RXD.COM")
            )
            typ = self.get_conds(rxd, "RXD.TYP")
            namedict = {
                ID: name
                for ID, name in zip(
                    coreinfo[0] + coreinfo[1] + rgt + cat + solv,
                    coreinfo[2] + coreinfo[3] + rgtnames + catnames + solvnames,
                )
            }
            yearpub = self.get_conds(rxd, "CIT.PREPY")
            nstages = self.get_conds(rxd, "RXD.SNR")
            yieldinfo = self.get_yield(rxd)
            if nstages:
                nstages = nstages[0]
            else:
                nstages = "1"
            rxndatadd = {
                "NumSteps": numsteps,
                "NumStages": nstages,
                "ReagentID": rgt,
                "MissingReagent": missingrgt,
                "Temperature": temp,
                "Pressure": press,
                "ReactionTime": time,
                "SolventID": solv,
                "MissingSolvent": missingsolv,
                "CatalystID": cat,
                "MissingCatalyst": missingcat,
                "ConditionNotes": condnotes,
                "ReactionType": typ,
                "NameDict": namedict,
                "YearPublished": yearpub,
                "Yield": yieldinfo,
            }
            details += [rxndatadd]

        # Do we expect for the current reaction more variations (RXD)?
        if len(rxd_list) < rxd_chunksize:
            return details

        #     breakpoint()
        # Exhaust the list of variations (RXD) for a reaction and create a result
        # as done for the initial slice above
        rxd_first = 1
        while True:
            rxd_first = rxd_first + rxd_chunksize
            rxd_last = rxd_first + rxd_chunksize - 1
            #         print(f"Processing RXDs from {rxd_first} to {rxd_last}")

            response = self.retrieve(
                self.resultname,
                [
                    "RX",
                    "RXD(%s,%s)" % (rxd_first, rxd_last),
                ],  # "RY","OMIT_MAPS=false,OMIT_CIT,HITONLY,ISSUE_RXN=false,ISSUE_RCT=false"
                current_rindex,
                current_rindex,
                "",
                "",
                "",
                "",
            )
            rxd_list = list(self.get_field_content(response, "RXD"))
            if len(rxd_list) == 0:
                break
            for rxd in rxd_list:
                numsteps = self.get_conds(rxd, "RXD.STP")[0]
                rgt, missingrgt, rgtnames = self.get_species(
                    rxd, "RXD.RGTXRN", "RXD.RGT", parenttag="RXD03"
                )
                cat, missingcat, catnames = self.get_species(
                    rxd, "RXD.CATXRN", "RXD.CAT", parenttag="RXD04"
                )
                solv, missingsolv, solvnames = self.get_species(
                    rxd, "RXD.SOLXRN", "RXD.SOL", parenttag="RXD05"
                )
                temp = self.get_conds(rxd, "RXD.T")
                press = self.get_conds(rxd, "RXD.P")
                time = self.get_conds(rxd, "RXD.TIM")
                condnotes = ", ".join(
                    self.get_conds(rxd, "RXD.COND") + self.get_conds(rxd, "RXD.COM")
                )
                typ = self.get_conds(rxd, "RXD.TYP")
                namedict = {
                    ID: name
                    for ID, name in zip(
                        coreinfo[0] + coreinfo[1] + rgt + cat + solv,
                        coreinfo[2] + coreinfo[3] + rgtnames + catnames + solvnames,
                    )
                }
                yearpub = self.get_conds(rxd, "CIT.PREPY")
                nstages = self.get_conds(rxd, "RXD.SNR")
                yieldinfo = self.get_yield(rxd)
                if nstages:
                    nstages = nstages[0]
                else:
                    nstages = "1"
                rxndatadd = {
                    "NumSteps": numsteps,
                    "NumStages": nstages,
                    "ReagentID": rgt,
                    "MissingReagent": missingrgt,
                    "Temperature": temp,
                    "Pressure": press,
                    "ReactionTime": time,
                    "SolventID": solv,
                    "MissingSolvent": missingsolv,
                    "CatalystID": cat,
                    "MissingCatalyst": missingcat,
                    "ConditionNotes": condnotes,
                    "ReactionType": typ,
                    "NameDict": namedict,
                    "YearPublished": yearpub,
                    "Yield": yieldinfo,
                }
                details += [rxndatadd]

        return details

    def get_reactions(
        self, hits_chunk, rxd_chunksize, rxids=[], rxid_first=None, rxid_last=None
    ):
        #     breakpoint()
        if rxids:
            self.select("RX", "R", f"RX.ID={';'.join(rxids)}", "", "WORKER,NO_CORESULT")
        elif rxid_first is not None and rxid_last is not None:
            self.select(
                "RX",
                "R",
                f"RXD.STP=1 AND RX.ID between {rxid_first} AND {rxid_last}",
                "",
                "WORKER,NO_CORESULT",
            )
        n_hits = int(self.resultsize)
        #     print(f"Resultsize: {n_hits}")
        rxndat = []
        # Keeps track of the currently processed reaction and the position within
        # the reaction hitset. This is required if additional reaction variations
        # have to be retrieved in case more than rxd_chunksize variations are available
        current_rindex = 0
        for rindex in range(1, n_hits + 1, hits_chunk):
            # Slice positions for reactions to retrieve
            first = rindex
            last = rindex + hits_chunk - 1
            #         print(f"From {first} to {last}")
            rxd_first = 1
            select_items = ["RX", "RXD(%s,%s)" % (rxd_first, rxd_chunksize)]  # "RY",
            response = self.retrieve(
                self.resultname, select_items, first, last, "", "", "", ""
            )  # "OMIT_MAPS=false,OMIT_CIT,HITONLY,ISSUE_RXN=false,ISSUE_RCT=false"
            #         breakpoint()
            reactions = self.get_field_content(response, "reaction")
            for reaction in reactions:
                # Keeping track of the current reaction position in the hitset
                current_rindex += 1
                rxndatcore = {}
                rxn_id = list(self.get_field_content(reaction, "RX.ID"))[
                    0
                ]  # Only one RX.ID per reaction available
                nref = list(self.get_field_content(reaction, "RX.NUMREF"))[0]
                rct, missingrct, rctnames = self.get_species(
                    reaction, "RX.RXRN", "RX.RCT", parenttag="RX01"
                )
                prod, missingprod, prodnames = self.get_species(
                    reaction, "RX.PXRN", "RX.PRO", parenttag="RX02"
                )
                rxndatcore.update(
                    {
                        "ReactionID": rxn_id,
                        "NumRefs": nref,
                        "ReactantID": rct,
                        "MissingReactant": missingrct,
                        "ProductID": prod,
                        "MissingProduct": missingprod,
                    }
                )
                #             breakpoint()
                details = self.get_details(
                    reaction,
                    [rct, prod, rctnames, prodnames],
                    current_rindex,
                    rxd_chunksize,
                )
                for detail in details:
                    #                 breakpoint()
                    variation = {**rxndatcore, **detail}
                    rxndat += [variation]
        return rxndat

    def main_(self, reaxys_ids=[], rxid_first=None, rxid_last=None):
        CALLERNAME = "chemeng_cam_201408"
        URL_MAIN = "www.reaxys.com"
        URL = "https://" + URL_MAIN + "/reaxys/api"
        # hit_chunks defines the number of reactions that will be returned in one
        # select/retrieve request. 200 is the recommended number
        hits_chunk = 100
        # A reaction can have multiple variations (RXD). To retrieve these variations
        # a slice must be defined in the retrieve statement as RXD(n,m). 20 is the
        # recommended number of RXDs to retrieve.
        rxd_chunksize = 20
        self.connect(URL, "", "", CALLERNAME)
        masterlst = []
        mastererrorlst = []
        if reaxys_ids:
            chunks = [
                reaxys_ids[i : i + hits_chunk]
                for i in range(0, len(reaxys_ids), hits_chunk)
            ]
            for chunk in chunks:
                try:
                    chunkres = self.get_reactions(
                        hits_chunk, rxd_chunksize, rxids=chunk
                    )
                    masterlst += chunkres
                except Exception as e:
                    mastererrorlst += chunk
        elif rxid_first is not None and rxid_last is not None:
            masterlst += self.get_reactions(
                hits_chunk, rxd_chunksize, rxid_first=rxid_first, rxid_last=rxid_last
            )
        self.disconnect()
        return masterlst, mastererrorlst


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import ray
from ray.util import ActorPool


def initray(restart=True, num_cpus=16, log_to_driver=False):
    """
    Initiates cluster of CPUs

    """
    if restart:
        ray.shutdown()
    ray.init(num_cpus=num_cpus, log_to_driver=log_to_driver)


def chunks(lst, s=None, k=None):
    """Yield successive s-sized chunks or k chunks from lst"""
    if s is not None:
        return [lst[i : i + s] for i in range(0, len(lst), s)]
    elif k is not None:
        n = len(lst)
        return [
            lst[i * (n // k) + min(i, n % k) : (i + 1) * (n // k) + min(i + 1, n % k)]
            for i in range(k)
        ]


def main__(num_cpus, reaxys_ids=[], rxid_first=None, rxid_last=None):
    initray(num_cpus=num_cpus)
    pool = ActorPool([Reaxys_API.remote() for i in range(num_cpus)])
    if reaxys_ids:
        masterchunks = chunks(
            reaxys_ids, k=num_cpus
        )  # Max 5 concurrent active sessions
        finalres = list(
            pool.map(lambda a, c: a.main_.remote(reaxys_ids=c), masterchunks)
        )
    elif rxid_first is not None and rxid_last is not None:
        rxid_first = min([rxid_first, rxid_last])
        rxid_last = max([rxid_first, rxid_last])
        masterchunks = chunks(
            list(range(int(rxid_first), int(rxid_last) + 1)), k=num_cpus
        )
        masterchunks = [[str(min(chunk)), str(max(chunk))] for chunk in masterchunks]
        finalres = list(
            pool.map(
                lambda a, c: a.main_.remote(rxid_first=c[0], rxid_last=c[1]),
                masterchunks,
            )
        )

    return finalres


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# rxn_ids=['639199'] #Specify list of IDs here. Note if reaction ID's to be retrieved by criteria, change where clause accordingly in get_reactions
# reaxys_final=[] #Final dictionary list
# while rxn_ids:
#     if len(rxn_ids)>=5:
#         nsessions=5 # Max 5 concurrent sessions
#     else:
#         nsessions=len(rxn_ids)
#     finalres=main__(nsessions,reaxys_ids=rxn_ids) #Calls API with 5 concurrent sessions, 100 IDs selected and retrieved at once
#     reaxys_dat=[record for batch in finalres for record in batch[0]]
#     errorlst=[errorid for batch in finalres for errorid in batch[1]]
#     if not reaxys_dat:
#         break
#     if not errorlst or len(errorlst)==1:
#         reaxys_final+=reaxys_dat
#         break
#     reaxys_final+=reaxys_dat
#     rxn_ids=errorlst
# # pd.to_pickle(reaxys_final,analoguedir+'analogue_rxns_updated_raw.pickle')
# # if errorlst:
# #     pd.to_pickle(errorlst,analoguedir+'error_list.pickle')
# #################Dataframe code######################
# analoguerxns_updated=pd.DataFrame(reaxys_final)
# for colname in ['ReactionID','NumRefs','NumSteps','NumStages']:
#     analoguerxns_updated[colname]=analoguerxns_updated[colname].apply(pd.to_numeric)
# for colname in ['ReactantID','ProductID','ReagentID','SolventID','CatalystID','YearPublished']:
#     analoguerxns_updated[colname]=analoguerxns_updated[colname].apply(lambda x: [int(ID) for ID in x])
# analoguerxns_updated['NameDict']=analoguerxns_updated['NameDict'].apply(lambda x: {int(key):val for key,val in x.items()})
# analoguerxns_updated['Yield']=analoguerxns_updated['Yield'].apply(lambda x: {int(key):float(val) for key,val in x.items()})
