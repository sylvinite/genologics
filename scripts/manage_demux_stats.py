#!/usr/bin/env python
"""
This file together with bcl_thresholds.py performs the bclconversion step of LIMS workflow.
In common tongue, it:
 
Fetches info from the sequencing process (RunID, FCID; derives instrument and data type)
Assigns (Q30, Clust per Lane) thresholds to the process (workflow step)
Reformats laneBarcode.html to "demuxstats_FCID_TIME.csv" for usage of other applications
Assigns a lot of info from laneBarcode.html to individual samples of the process (%PF etc)
Flags samples as QC PASS/FAIL based on thresholds

##As of implementation preproc2 has hiseq examples & Preproc 1 has hiseqX
##lanebc file is stored in Demultiplexing_0 on nosync on /srv/illumina

Written by Isak Sylvin; isak.sylvin@scilifelab.se"""

from genologics.lims import Lims
from genologics.config import BASEURI, USERNAME, PASSWORD
from genologics.entities import Process
from genologics.epp import EppLogger, attach_file
import flowcell_parser.classes as classes
from bcl_thresholds import Thresholds
from datetime import datetime
from time import time

import os 
import click
import csv
import sys
import logging

timestamp = datetime.fromtimestamp(time()).strftime("%Y-%m-%d_%H:%M")

"""Fetches overarching workflow info"""
def manipulate_workflow(demux_process):
    try:
        demux_container = demux_process.all_inputs()[0].location[0].name    
    except:
        sys.exit("Container name not found")
        
    workflow_types = {"MiSeq Run (MiSeq) 4.0":"Reagent Cartridge ID", "Illumina Sequencing (Illumina SBS) 4.0":"Flow Cell ID",
              "Illumina Sequencing (HiSeq X) 1.0":"Flow Cell ID"}
    for k,v in workflow_types.items():
        #Sequencing step is null if it doesn"t exist
        workflow = lims.get_processes(udf = {v : demux_container}, type = k)
        #If sequencing step key exists
        if(workflow):
            #Copies LIMS sequencing step content
            proc_stats = dict(workflow[0].udf.items())
            #Instrument is denoted the way it is since it is also used to find
            #the folder of the laneBarcode.html file
            if "MiSeq Run (MiSeq) 4.0" in k:
                proc_stats["Chemistry"] ="MiSeq"
                proc_stats["Instrument"] = "miseq"
            elif "Illumina Sequencing (Illumina SBS) 4.0" in k:
                try:
                    proc_stats["Chemistry"] = workflow[0].udf["Flow Cell Version"]
                except:
                    sys.exit("No flowcell version set in sequencing step.")
                proc_stats["Instrument"] = "hiseq"
            elif "Illumina Sequencing (HiSeq X) 1.0" in k:
                proc_stats["Chemistry"] ="HiSeqX v2.5"
                proc_stats["Instrument"] = "hiseqx"
            else:
                sys.exit("Unhandled prior workflow step (run type)")
            logging.info("Run type set to {}".format(proc_stats["Chemistry"]))
            break
    
    try:
        proc_stats["Paired"] = False
    except:
        sys.exit("Unable to fetch workflow information.")
    if "Read 2 Cycles" in proc_stats:
        proc_stats["Paired"] = True
    logging.info("Paired libraries: {}".format(str(proc_stats["Paired"])))  
    #Assignment to make usage more explicit
    try:
        proc_stats["Read Length"] = proc_stats["Read 1 Cycles"]
    except:
        sys.exit("Read 1 Cycles not found. Unable to read Read Length")
    return proc_stats

"""Sets run thresholds"""
def manipulate_process(demux_process, proc_stats):      
    thresholds = Thresholds(proc_stats["Instrument"], proc_stats["Chemistry"], proc_stats["Paired"], proc_stats["Read Length"])

    if not "Threshold for % bases >= Q30" in demux_process.udf:
        try:
            demux_process.udf["Threshold for % bases >= Q30"] = thresholds.Q30
        except:
            sys.exit("Udf improperly formatted. Unable to set Q30 threshold")
    #Would prefer "Expected reads per lane" but "Threshold for # Reads" means less LIMS alterations
    if not "Threshold for # Reads" in demux_process.udf:
        try:
            demux_process.udf["Threshold for # Reads"] = thresholds.exp_lane_clust
        except:
            sys.exit("Udf improperly formatted. Unable to set # Reads threshold")
    
    #Would REALLY prefer "Threshold for % undetermined reads per lane"
    if not "Threshold for Undemultiplexed Index Yield" in demux_process.udf:
        try:
            demux_process.udf["Threshold for Undemultiplexed Index Yield"] = thresholds.undet_indexes_perc
        except:
            sys.exit("Udf improperly formatted. Unable to set Undemultiplexed Index Yield threshold")

    logging.info("Q30 threshold set to {}".format(str(demux_process.udf["Threshold for % bases >= Q30"])))
    logging.info("Minimum clusters per lane set to {}".format(str(demux_process.udf["Threshold for # Reads"])))
    logging.info("Maximum percentage of undetermined per lane set to {} %".format(str(demux_process.udf["Threshold for Undemultiplexed Index Yield"])))
    
    
    #Sets run id if not already exists:
    if not "Run ID" in demux_process.udf:
        try:
            demux_process.udf["Run ID"] = proc_stats["Run ID"]
        except:
            logging.info("Unable to automatically regenerate Run ID")
    #Check for document version
    if not "Document Version" in demux_process.udf:
        sys.exit("No Document Version set. Please set one.")
    try:
        demux_process.put()
    except:
        sys.exit("Failed to apply process thresholds to LIMS")
    
"""Sets artifact = samples values """
def set_sample_values(demux_process, parser_struct, proc_stats):
    for pool in demux_process.all_inputs():
        try:
            outarts_per_lane = demux_process.outputs_per_input(pool.id, ResultFile = True)
        except:
            sys.exit("Unable to fetch artifacts of process")
        if proc_stats["Instrument"] == "miseq":
            lane_no = "1"
        else:
            try:
                lane_no = pool.location[1][0]
            except:
                sys.exit("Unable to determine lane number. Incorrect location variable in process.")
        logging.info("Lane number set to {}".format(lane_no))
        exp_smp_per_lne = round(demux_process.udf["Threshold for # Reads"]/float(len(outarts_per_lane)), 0)
        logging.info("Expected sample clusters for this lane: {}".format(str(exp_smp_per_lne)))
        for target_file in outarts_per_lane:
            try:
                current_name = target_file.samples[0].name
            except:
                sys.exit("Unable to determine sample name. Incorrect sample variable in process.")
            for entry in parser_struct:
                if lane_no == entry["Lane"]:
                    total_reads = 0
                    sample = entry["Sample"]
                    #TODO: Add lane warning
                    if sample == current_name:
                        try:
                            target_file.udf["%PF"] = float(entry["% PFClusters"])
                            target_file.udf["% One Mismatch Reads (Index)"] = float(entry["% One mismatchbarcode"])
                            target_file.udf["% of Raw Clusters Per Lane"] = float(entry["% of thelane"])
                            target_file.udf["Ave Q Score"] = float(entry["Mean QualityScore"])
                            target_file.udf["% Perfect Index Read"] = float(entry["% Perfectbarcode"])
                            target_file.udf["Yield PF (Gb)"] = float(entry["Yield (Mbases)"].replace(",",""))/1000
                            target_file.udf["% Bases >=Q30"] = float(entry["% >= Q30bases"])
                            #Paired runs are divided by two within flowcell parser
                            if proc_stats["Paired"]:
                                target_file.udf["# Reads"] = int(entry["PF Clusters"].replace(",",""))*2
                                target_file.udf["# Read Pairs"] = int(entry["PF Clusters"].replace(",",""))
                            #Since a single ended run has no pairs, pairs is set to equal reads
                            else:
                                target_file.udf["# Reads"] = int(entry["PF Clusters"].replace(",",""))
                                target_file.udf["# Read Pairs"] = int(entry["PF Clusters"].replace(",",""))
                            total_reads = total_reads + target_file.udf["# Reads"]
                        except:
                            sys.exit("Unable to set artifact values.")
                        logging.info("Added the following set of values to {} of lane {}:".format(sample,lane_no))
                        logging.info("{}% PF".format(str(target_file.udf["%PF"])))
                        logging.info("{}% One Mismatch Reads (Index)".format(str(target_file.udf["% One Mismatch Reads (Index)"])))
                        logging.info("{}% of Raw Clusters Per Lane".format(str(target_file.udf["% of Raw Clusters Per Lane"])))
                        logging.info("{}Ave Q Score".format(str(target_file.udf["Ave Q Score"])))
                        logging.info("{}% Perfect Index Read".format(str(target_file.udf["% Perfect Index Read"])))
                        logging.info("{} Yield (Mbases)".format(str(target_file.udf["Yield PF (Gb)"])))
                        logging.info("{}% Bases >=Q30".format(str(target_file.udf["% Bases >=Q30"])))
                        logging.info("{}# Reads".format(str(target_file.udf["# Reads"])))
                        logging.info("{}# Reads Pairs".format(str(target_file.udf["# Read Pairs"])))
                        
                        try:
                            if (demux_process.udf["Threshold for % bases >= Q30"] <= float(entry["% >= Q30bases"]) and 
                            int(exp_smp_per_lne) <= target_file.udf["# Reads"] ):
                                target_file.udf["Include reads"] = "YES"
                                target_file.qc_flag = "PASSED"           
                            else:
                                target_file.udf["Include reads"] = "NO"
                                target_file.qc_flag = "FAILED"
                            logging.info("Q30 %: {} versus {}".format(str(float(entry["% >= Q30bases"])), str(demux_process.udf["Threshold for % bases >= Q30"])))
                            logging.info("Expected reads: {} versus {}".format(str(target_file.udf["# Reads"]), str(int(exp_smp_per_lne)))) 
                            logging.info("Sample QC status set to {}".format(target_file.qc_flag))
                        except:
                            sys.exit("Unable to set QC status for sample")
                #TODO: REVIEW THIS WHOLE SECTION
                    elif sample == "Undetermined":
                        if proc_stats["Paired"]:
                                reads = int(entry["PF Clusters"].replace(",",""))*2
                                read_pairs= int(entry["PF Clusters"].replace(",",""))
                        else:
                                reads = int(entry["PF Clusters"].replace(",",""))
                                read_pairs = int(entry["PF Clusters"].replace(",",""))
                        logging.info("Undemultiplexed reads at {} % for lane {}".format(reads/(total_reads + reads), lane_no))
                        #If undetermined reads are greater than %*reads_in_lane
                        if reads/(total_reads + reads) > demux_process.udf["Threshold for Undemultiplexed Index Yield"]:
                            sys.stderr("WARNING: Undemultiplexed reads for lane {} exceeded defined limit".format(lane_no))
            try: 
                target_file.put()
            except:
                sys.exit("Failed to apply artifact data to LIMS")
    

"""Creates demux_{FCID}_{time}.csv and attaches it to process"""
def write_demuxfile(proc_stats):
    #Includes windows drive letter support
    datafolder = "{}_data".format(proc_stats["Instrument"])
    lanebc_path = os.path.join(os.sep,"srv","mfs", "taca_test",datafolder,proc_stats["Run ID"],"laneBarcode.html")
    try:
        laneBC = classes.LaneBarcodeParser(lanebc_path)
    except:
        sys.exit("Unable to fetch laneBarcode.html from {}".format(lanebc_path))
    fname = "demuxstats_{}_{}.csv".format(proc_stats["Flow Cell ID"], timestamp)
    
    #Dumps less undetermined info than undemultiplex_index.py. May cause problems downstreams
    with open(fname, "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Project", "Sample ID", "Lane", "# Reads", "Index","Index name", "% of >= Q30 Bases (PF)"])
        for entry in laneBC.sample_data:
            index_name = ""
            if "Clusters" in entry:
                reads = entry["Clusters"]
            else: 
                reads = entry["PF Clusters"]
                
            if proc_stats["Paired"]:
                reads = int(reads.replace(",",""))*2
            else:
                reads = int(reads.replace(",","")) 
            
            try:
                writer.writerow([entry["Project"],entry["Sample"],entry["Lane"],reads,entry["Barcode sequence"],index_name,entry["% >= Q30bases"]])
            except:
                sys.exit("Flowcell parser is unable to fetch all necessary fields for demux file.")
    return laneBC.sample_data

#Debug command (p1)
#def converter(demux_process):
def converter(demux_process, epp_logger):

    #Fetches info on "workflow" level
    proc_stats = manipulate_workflow(demux_process)
    #Sets up the process values
    manipulate_process(demux_process, proc_stats)
    #Create the demux output file
    parser_struct = write_demuxfile(proc_stats)
    #Alters artifacts
    set_sample_values(demux_process, parser_struct, proc_stats)
    
    #Attaches output files to lims process; LIMS magic
    for out in demux_process.all_outputs():
        if out.name == "Demultiplex Stats":
            attach_file("demuxstats_{}_{}.csv".format(proc_stats["Flow Cell ID"], timestamp), out)
        elif out.name == "Script Log Details":
            attach_file(epp_logger.log_file, out)

@click.command()
@click.option("--project_lims_id", required=True,help="REQUIRED: Lims ID of project. Example:24-92373")
@click.option("--rt_log", help="File name to log runtime information.") 
def main(project_lims_id, rt_log):
    demux_process = Process(lims,id = project_lims_id)
    #Debug command (p2)
    #converter(demux_process)
    #Sets up proper logging
    with EppLogger(log_file=rt_log, lims=lims, prepend=True) as epp_logger:
        converter(demux_process, epp_logger)      
if __name__ == "__main__":
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main()
