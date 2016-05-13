#!/usr/bin/env python
"""
This file together with manage_demux_stats_thresholds.py performs the "bclconversion" step of LIMS workflow.
In common tongue, it:
 
Fetches info from the sequencing process (RunID, FCID; derives instrument and data type)
Assigns (Q30, Clust per Lane) thresholds to the process (workflow step)
Reformats laneBarcode.html to "demuxstats_FCID.csv" for usage of other applications
Assigns a lot of info from laneBarcode.html to individual samples of the process (e.g. %PF)
Flags samples as QC PASSED/FAILED based on thresholds

Written by Isak Sylvin; isak.sylvin@scilifelab.se"""

#Fetched from SciLife repos
from genologics.lims import Lims
from genologics.config import BASEURI, USERNAME, PASSWORD
from genologics.entities import Process
from genologics.epp import EppLogger, attach_file
import flowcell_parser.classes as classes
from manage_demux_stats_thresholds import Thresholds

#Standard packages
from shutil import move
import os 
import click
import csv
import sys
import logging

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
            logging.info("Q30 threshold set to {}".format(str(demux_process.udf["Threshold for % bases >= Q30"])))
        except:
            sys.exit("Udf improperly formatted. Unable to set Q30 threshold")
    #Would REALLY prefer "Minimum Reads per lane" over "Threshold for # Reads"
    if not "Minimum Reads per lane" in demux_process.udf:
        try:
            demux_process.udf["Minimum Reads per lane"] = thresholds.exp_lane_clust
            logging.info("Minimum clusters per lane set to {}".format(str(demux_process.udf["Minimum Reads per lane"])))
        except:
            sys.exit("Udf improperly formatted. Unable to set # Reads threshold")
    
    #Would REALLY prefer "Maximum % Undetermined Reads per Lane" over "Threshold for Undemultiplexed Index Yield"
    if not "Maximum % Undetermined Reads per Lane" in demux_process.udf:
        try:
            demux_process.udf["Maximum % Undetermined Reads per Lane"] = thresholds.undet_indexes_perc
            logging.info("Maximum percentage of undetermined per lane set to {} %".format(str(demux_process.udf["Maximum % Undetermined Reads per Lane"])))
        except:
            sys.exit("Udf improperly formatted. Unable to set Undemultiplexed Index Yield threshold")

    #Sets Run ID if not already exists:
    if not "Run ID" in demux_process.udf:
        try:
            demux_process.udf["Run ID"] = proc_stats["Run ID"]
        except:
            logging.info("Unable to automatically regenerate Run ID")
    #Checks for document version
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
        exp_smp_per_lne = round(demux_process.udf["Minimum Reads per lane"]/float(len(outarts_per_lane)), 0)
        logging.info("Expected sample clusters for this lane: {}".format(str(exp_smp_per_lne)))
        
        assign_lane_reads = 0
        undet_lane_reads = 0
        for target_file in outarts_per_lane:
            try:
                current_name = target_file.samples[0].name
            except:
                sys.exit("Unable to determine sample name. Incorrect sample variable in process.")
            for entry in parser_struct:
                if lane_no == entry["Lane"]:
                    sample = entry["Sample"]
                    if sample == current_name:
                        logging.info("Added the following set of values to {} of lane {}:".format(sample,lane_no))
                        try:
                            target_file.udf["%PF"] = float(entry["% PFClusters"])
                            logging.info("{}% PF".format(str(target_file.udf["%PF"])))
                            target_file.udf["% One Mismatch Reads (Index)"] = float(entry["% One mismatchbarcode"])
                            logging.info("{}% One Mismatch Reads (Index)".format(str(target_file.udf["% One Mismatch Reads (Index)"])))
                            target_file.udf["% of Raw Clusters Per Lane"] = float(entry["% of thelane"])
                            logging.info("{}% of Raw Clusters Per Lane".format(str(target_file.udf["% of Raw Clusters Per Lane"])))
                            target_file.udf["Ave Q Score"] = float(entry["Mean QualityScore"])
                            logging.info("{}Ave Q Score".format(str(target_file.udf["Ave Q Score"])))
                            target_file.udf["% Perfect Index Read"] = float(entry["% Perfectbarcode"])
                            logging.info("{}% Perfect Index Read".format(str(target_file.udf["% Perfect Index Read"])))
                            target_file.udf["Yield PF (Gb)"] = float(entry["Yield (Mbases)"].replace(",",""))/1000
                            logging.info("{} Yield (Mbases)".format(str(target_file.udf["Yield PF (Gb)"])))
                            target_file.udf["% Bases >=Q30"] = float(entry["% >= Q30bases"])
                            logging.info("{}% Bases >=Q30".format(str(target_file.udf["% Bases >=Q30"])))
                        except:
                            sys.exit("Unable to set general artifact values")
                        try:
                            clusterType = None
                            if "PF Clusters" in entry:
                                clusterType = "PF Clusters"
                            else:
                                clusterType = "Clusters"
                            #Paired runs are divided by two within flowcell parser
                            if proc_stats["Paired"]:
                                target_file.udf["# Reads"] = int(entry[clusterType].replace(",",""))*2
                                target_file.udf["# Read Pairs"] = int(entry[clusterType].replace(",",""))
                            #Since a single ended run has no pairs, pairs is set to equal reads
                            else:
                                target_file.udf["# Reads"] = int(entry[clusterType].replace(",",""))
                                target_file.udf["# Read Pairs"] = int(entry[clusterType].replace(",",""))
                            assign_lane_reads = assign_lane_reads + target_file.udf["# Reads"]
                        except:
                            sys.exit("Unable to set values for #Reads and #Read Pairs.")
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
                            logging.info("Q30 %: {}% found, minimum at {}%".format(str(float(entry["% >= Q30bases"])), str(demux_process.udf["Threshold for % bases >= Q30"])))
                            logging.info("Expected reads: {} found, minimum at {}".format(str(target_file.udf["# Reads"]), str(int(exp_smp_per_lne)))) 
                            logging.info("Sample QC status set to {}".format(target_file.qc_flag))
                        except:
                            sys.exit("Unable to set QC status for sample")
                    #Counts undetermined per lane
                    elif sample == "Undetermined":
                        clusterType = None
                        if "PF Clusters" in entry:
                            clusterType = "PF Clusters"
                        else:
                            clusterType = "Clusters"
                            
                        if proc_stats["Paired"]:
                            undet_lane_reads = int(entry[clusterType].replace(",",""))*2
                            undet_read_pairs= int(entry[clusterType].replace(",",""))
                        else:
                            undet_lane_reads = int(entry[clusterType].replace(",",""))
                            undet_read_pairs = int(entry[clusterType].replace(",",""))
                        
            try: 
                target_file.put()
            except:
                sys.exit("Failed to apply artifact data to LIMS")
            
        #If undetermined reads are greater than threshold*reads_in_lane
        found_undet = round(float(undet_lane_reads)/(assign_lane_reads+undet_lane_reads)*100, 2)
        if  found_undet > demux_process.udf["Maximum % Undetermined Reads per Lane"]:
            #This looks kind of bad in lims, but \n can't be used.
            sys.stderr.write("WARNING: Undemultiplexed reads for lane {} was {} ({})% thus exceeding defined limit.\t".format(lane_no, undet_lane_reads, found_undet))
        else:
            logging.info("Found {} ({}%) undemultiplexed reads for lane {}.".format(undet_lane_reads, found_undet, lane_no))

"""Creates demux_{FCID}.csv and attaches it to process"""
def write_demuxfile(proc_stats, demux_id):
    #Includes windows drive letter support
    datafolder = "{}_data".format(proc_stats["Instrument"])
    lanebc_path = os.path.join(os.sep,"srv","mfs", "taca_test",datafolder,proc_stats["Run ID"],"laneBarcode.html")
    try:
        laneBC = classes.LaneBarcodeParser(lanebc_path)
    except:
        sys.exit("Unable to fetch laneBarcode.html from {}".format(lanebc_path))
    fname = "{}_demuxstats_{}.csv".format(demux_id, proc_stats["Flow Cell ID"])
    
    #Writes less undetermined info than undemultiplex_index.py. May cause problems downstreams
    with open(fname, "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Project", "Sample ID", "Lane", "# Reads", "Index","Index name", "% of >= Q30 Bases (PF)"])
        for entry in laneBC.sample_data:
            index_name = ""
            if "PF Clusters" in entry:
                reads = entry["PF Clusters"]
            else: 
                reads = entry["Clusters"]
                
            if proc_stats["Paired"]:
                reads = int(reads.replace(",",""))*2
            else:
                reads = int(reads.replace(",","")) 
            
            try:
                writer.writerow([entry["Project"],entry["Sample"],entry["Lane"],reads,entry["Barcode sequence"],index_name,entry["% >= Q30bases"]])
            except:
                sys.exit("Flowcell parser is unable to fetch all necessary fields for demux file.")
    return laneBC.sample_data

def converter(demux_process, epp_logger,demux_id):
    
    #Fetches info on "workflow" level
    proc_stats = manipulate_workflow(demux_process)
    #Sets up the process values
    manipulate_process(demux_process, proc_stats)
    #Create the demux output file
    parser_struct = write_demuxfile(proc_stats, demux_id)
    #Alters artifacts
    set_sample_values(demux_process, parser_struct, proc_stats)
    
    #Changing log file name, can't do this step earlier since it breaks logging
    for out in demux_process.all_outputs():
        if out.name == "Script Log Details":
            new_name = "{}_logfile_{}.txt".format(epp_logger.log_file, proc_stats["Flow Cell ID"])
            move(epp_logger.log_file, new_name)

@click.command()
@click.option("--project_lims_id", required=True,help="Lims ID of project. Example:24-92373")
@click.option("--demux_id", required=True,help="Id prefix for demux output.")
@click.option("--log_id", required=True,help="Id prefix for logfile")
def main(project_lims_id, demux_id, log_id ):
    demux_process = Process(lims,id = project_lims_id)
    #Sets up proper logging
    with EppLogger(log_file=log_id, lims=lims, prepend=True) as epp_logger:
        converter(demux_process, epp_logger, demux_id)      
if __name__ == "__main__":
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main()
