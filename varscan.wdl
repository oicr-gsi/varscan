version 1.0

struct GenomeResources {
    String workflowModules
    String refFasta
    String mergingModules
    String refDict
    Int largest
}

workflow varscan {
input {
    # Normally we need only tumor bam, normal bam may be used when available
    File inputTumor
    File inputNormal
    File inputTumorIndex
    File inputNormalIndex
    String reference
    String outputFileNamePrefix
    String bedIntervalsPath = ""
}

Map[String, Array[String]] regions = {
  "hg19": ["chr1:1-249250621","chr2:1-243199373","chr3:1-198022430","chr4:1-191154276","chr5:1-180915260","chr6:1-171115067","chr7:1-159138663","chr8:1-146364022","chr9:1-141213431","chr10:1-135534747","chr11:1-135006516","chr12:1-133851895","chr13:1-115169878","chr14:1-107349540","chr15:1-102531392","chr16:1-90354753","chr17:1-81195210","chr18:1-78077248","chr19:1-59128983","chr20:1-63025520","chr21:1-48129895","chr22:1-51304566","chrX:1-155270560","chrY:1-59373566","chrM:1-16571"],
  "hg38": ["chr1:1-248956422","chr2:1-242193529","chr3:1-198295559","chr4:1-190214555","chr5:1-181538259","chr6:1-170805979","chr7:1-159345973","chr8:1-145138636","chr9:1-138394717","chr10:1-133797422","chr11:1-135086622","chr12:1-133275309","chr13:1-114364328","chr14:1-107043718","chr15:1-101991189","chr16:1-90338345","chr17:1-83257441","chr18:1-80373285","chr19:1-58617616","chr20:1-64444167","chr21:1-46709983","chr22:1-50818468","chrX:1-156040895","chrY:1-57227415","chrM:1-16569"],
  "mm10": ["chr1:1-195471971","chr2:1-182113224","chr3:1-160039680","chr4:1-156508116","chr5:1-151834684","chr6:1-149736546","chr7:1-145441459","chr8:1-129401213","chr9:1-124595110","chr10:1-130694993","chr11:1-122082543","chr12:1-120129022","chr13:1-120421639","chr14:1-124902244","chr15:1-104043685","chr16:1-98207768","chr17:1-94987271","chr18:1-90702639","chr19:1-61431566","chrX:1-171031299","chrY:1-91744698","chrM:1-16299"]
}

Map[String,GenomeResources] resources = {
  "hg19":  {
    "workflowModules": "samtools/0.1.19 hg19/p13",
    "refFasta": "$HG19_ROOT/hg19_random.fa",
    "mergingModules": "picard/2.21.2 hg19/p13",
    "refDict": "$HG19_ROOT/hg19_random.dict",
    "largest": "249250621"
  },
  "hg38": {
    "workflowModules": "samtools/0.1.19 hg38/p12",
    "refFasta": "$HG38_ROOT/hg38_random.fa",
    "mergingModules": "picard/2.21.2 hg38/p12",
    "refDict": "$HG38_ROOT/hg38_random.dict",
    "largest": "248956422"
  },
  "mm10": {
    "workflowModules": "samtools/0.1.19 mm10/p6",
    "refFasta": "$MM10_ROOT/mm10.fa",
    "mergingModules": "picard/2.21.2 mm10/p6",
    "refDict": "$MM10_ROOT/mm10.dict",
    "largest": "195471971"
  } 
}


call expandRegions { input: bedPath = bedIntervalsPath }

String sampleID = if outputFileNamePrefix=="" then basename(inputTumor, ".bam") else outputFileNamePrefix
Array[String] splitRegions = if bedIntervalsPath != "" then expandRegions.regions else regions[reference]

# Produce pileups
scatter ( r in splitRegions )   {
  call getChrCoefficient {
      input: modules = resources [ reference ].workflowModules,
             refDict = resources [ reference ].refDict,
             largestChrom = resources [ reference ].largest,
             region = r
  }

  call makePileups { input: inputTumor = inputTumor,
                            inputTumorIndex = inputTumorIndex,
                            inputNormal = inputNormal,
                            inputNormalIndex = inputNormalIndex,
                            modules = resources[reference].workflowModules,
                            refFasta = resources[reference].refFasta,
                            scaleCoefficient = getChrCoefficient.coeff,
                            region = r }

  # Configure and run Varscan
  call runVarscanCNV { input: inputPileup = makePileups.pileup, sampleID = sampleID, scaleCoefficient = getChrCoefficient.coeff }
  call runVarscanSNV as getSnvNative { input: inputPileup = makePileups.pileup, sampleID = sampleID, scaleCoefficient = getChrCoefficient.coeff }
  call runVarscanSNV as getSnvVcf { input: inputPileup = makePileups.pileup, sampleID = sampleID, outputVcf = 1, scaleCoefficient = getChrCoefficient.coeff }
}

# Merge tasks
call mergeVariantsNative as mergeCNV { input: filePaths = select_all(runVarscanCNV.resultFile), outputFile = sampleID, outputExtension = "copynumber" }
call mergeVariantsNative as mergeSNP { input: filePaths = select_all(getSnvNative.snpFile), outputFile = sampleID, outputExtension = "snp" }
call mergeVariantsNative as mergeIND { input: filePaths = select_all(getSnvNative.indelFile), outputFile = sampleID, outputExtension = "indel" }
call mergeVariantsVcf as mergeSNPvcf { input: filePaths = select_all(getSnvVcf.snpVcfFile), 
                                              outputSuffix = "snp",
                                              outputFile = sampleID,
                                              seqDictionary = resources[reference].refDict,
                                              modules = resources[reference].mergingModules }
call mergeVariantsVcf as mergeINDvcf { input: filePaths = select_all(getSnvVcf.indelVcfFile),
                                              outputSuffix = "indel",
                                              outputFile = sampleID,
                                              seqDictionary = resources[reference].refDict,
                                              modules = resources[reference].mergingModules }

# Run post-processing job if we have results from runVarscanCNV
Array[File] cNumberFile = [mergeCNV.mergedVariants]
if (length(cNumberFile) == 1) {
    call smoothData{input: copyNumberFile = mergeCNV.mergedVariants, sampleID = sampleID}
}

### merge the snv and indel vcf into a single file
call vcfCombine {   input: vcfSnvs = mergeSNPvcf.mergedVcf,
                           vcfIndels = mergeINDvcf.mergedVcf,
                           outputFileNamePrefix = outputFileNamePrefix
}

meta {
  author: "Peter Ruzanov"
  email: "pruzanov@oicr.on.ca"
  description: "Varscan 2.3, workflow for calling SNVs and CVs\nCreation of mpileups and calling variants are done with parallel processing\n\n![varscan outputs](docs/Screenshot_Varscan.png)"

  dependencies: [
      {
        name: "picard/2.21.2",
        url: "https://broadinstitute.github.io/picard"
      },
      {
        name: "varscan/2.4.2",
        url: "http://varscan.sourceforge.net"
      },
      {
        name: "samtools/0.1.19",
        url: "http://www.htslib.org/"
      },
      {
        name: "rstats/3.6",
        url: "http://cran.utstat.utoronto.ca/src/base/R-3/R-3.6.1.tar.gz"
      }
    ]
    
    output_meta: {
    resultCnvFile: {
        description: "file with CNV calls, smoothed",
        vidarr_label: "resultCnvFile"
    },
    resultSnpFile: {
        description: "file with SNPs, native varscan format",
        vidarr_label: "resultSnpFile"
    },
    resultIndelFile: {
        description: "file with Indel calls, native varscan format",
        vidarr_label: "resultIndelFile"
    },
    resultSnpVcfFile: {
        description: "file with SNPs, vcf format",
        vidarr_label: "resultSnpVcfFile"
    },
    resultIndelVcfFile: {
        description: "file with Indels, vcf format",
        vidarr_label: "resultIndelVcfFile"
    },
    resultVcfFile: {
        description: "file with snvs + indels, vcf format, bgzipped",
        vidarr_label: "resultVcfFile"
    },
    resultVcfFileIndex: {
        description: "index file for snv + indels vcf output",
        vidarr_label: "resultVcfFileIndex"
    }


}
}

parameter_meta {
  inputTumor: "input .bam file for tumor sample"
  inputNormal: "input .bam file for normal sample"
  inputTumorIndex: "input .bai file for tumor sample"
  inputNormalIndex: "input .bai file for normal sample"
  reference: "Reference assembly id, hg19 hg38 or mm10"
  outputFileNamePrefix: "Output file(s) prefix"
  bedIntervalsPath: "Path to a .bed file used for splitting pileup job/limiting analysis to selected regions"
  chromRegions: "Regions used for scattering tasks, need to be assembly-specific"
}

output {
 File? resultCnvFile     = smoothData.filteredData
 File resultSnpFile      = mergeSNP.mergedVariants
 File resultIndelFile    = mergeIND.mergedVariants
 File resultSnpVcfFile   = mergeSNPvcf.mergedVcf
 File resultIndelVcfFile = mergeINDvcf.mergedVcf
 File resultVcfFile      = vcfCombine.vcf
 File resultVcfFileIndex = vcfCombine.vcfIndex
}

}


# ================================================================
#  Combine the indel and snv vcfs into a single ouutput
# ================================================================

task vcfCombine {
    input {
	   File vcfSnvs
	   File vcfIndels
	   String modules = "bcftools/1.9 tabix/1.9"
	   String outputFileNamePrefix
	   Int jobMemory = 16
	   Int threads = 4
	   Int timeout = 4	   
	
	}

    parameter_meta {
	vcfSnvs: "vcf file with snvs"
	vcfIndels: "vcf file with indels"
	modules: "environment modules"
	jobMemory: "Memory allocated for job"
	timeout: "Hours before task timeout"
	threads: "Number of threads for processing"
    }

	
    command <<<
	set -eo pipefail
	
	bgzip -c ~{vcfSnvs} > ~{outputFileNamePrefix}.varscan2_snv.vcf.gz
	tabix ~{outputFileNamePrefix}.varscan2_snv.vcf.gz
	bgzip -c ~{vcfIndels} > ~{outputFileNamePrefix}.varscan2_indel.vcf.gz
	tabix ~{outputFileNamePrefix}.varscan2_indel.vcf.gz
	bcftools concat -a -o ~{outputFileNamePrefix}.varscan2_all.vcf ~{outputFileNamePrefix}.varscan2_snv.vcf.gz ~{outputFileNamePrefix}.varscan2_indel.vcf.gz
	bgzip ~{outputFileNamePrefix}.varscan2_all.vcf
	tabix ~{outputFileNamePrefix}.varscan2_all.vcf.gz

	>>>
	
    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File vcf = "~{outputFileNamePrefix}.varscan2_all.vcf.gz"
	File vcfIndex = "~{outputFileNamePrefix}.varscan2_all.vcf.gz.tbi"
    }
    
	meta {
	output_meta: {
	    vcf: "VCF file with snvs and indels, bgzip compressed",
	    vcfIndex: "tabix index"
	}
    }
		

}





# ================================================================
#  Scaling coefficient - use to scale RAM allocation by chromosome
# ================================================================
task getChrCoefficient {
  input {
    Int memory = 1
    Int timeout = 1
    Int largestChrom
    String region
    String modules
    String refDict
  }

  parameter_meta {
    refDict: ".dict file for the reference genome, we use it to extract chromosome ids"
    timeout: "Hours before task timeout"
    region: "Region to extract a chromosome to check"
    memory: "Memory allocated for this job"
    modules: "Names and versions of modules to load"
    largestChrom: "Length of the largest chromosome in a genome"
  }

  command <<<
    CHROM=$(echo ~{region} | sed 's/:.*//')
    grep -w SN:$CHROM ~{refDict} | cut -f 3 | sed 's/.*://' | awk '{print int(($1/~{largestChrom} + 0.1) * 10)/10}'
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    String coeff = read_string(stdout())
  }

  meta {
    output_meta: {
      coeff: "Length ratio as relative to the largest chromosome."
    }
  }
}

# =======================================================
# Read bed file, return a string with regions for mpileup
# =======================================================
task expandRegions {
input {
 String bedPath = ""
 String modules = "hg38-dac-exclusion/1.0"
 Int jobMemory = 4
 Int timeout = 12
}

parameter_meta {
  bedPath: "Optional path to a bed file with intervals for splitting pileup job/limiting to regions"
  jobMemory: "Memory for this task in GB"
  modules: "required modules (This is to allow modularized data for bed path)" 
  timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
 python <<CODE
 import os
 if os.path.exists("~{bedPath}"):
    with open("~{bedPath}") as f:
        for line in f:
           line = line.rstrip()
           tmp = line.split("\t")
           r = " " + tmp[0] + ":" + tmp[1] + "-" + tmp[2]
           print(r)
    f.close()
 CODE
>>>

runtime {
 memory:  "~{jobMemory} GB"
 modules: "~{modules}"
 timeout: "~{timeout}"
}

output {
 Array[String] regions = read_lines(stdout()) 
}
}

# ==========================================
#  produce pileup with samtools
# ==========================================
task makePileups {
input {
 File inputNormal
 File inputTumor
 File inputTumorIndex
 File inputNormalIndex
 String refFasta 
 String modules
 String samtools = "$SAMTOOLS_ROOT/bin/samtools"
 String region 
 Float scaleCoefficient = 1.0
 Int jobMemory   = 18
 Int minMemory = 6
 Int timeout     = 40
}

parameter_meta {
  inputNormal: "input .bam file for normal tissue"
  inputNormalIndex: ".bai index file for normal tissue"
  inputTumor: "input .bam file for tumor tissue"
  inputTumorIndex: ".bai index file for tumor tissue"
  refFasta: "Reference fasta file, path depends on the respective module"
  modules: "required modules"
  samtools: "path to samtools"
  region: "Region in a form of chrX:12000-12500 for mpileup command"
  scaleCoefficient: "Scaling coefficient for RAM allocation, depends on chromosome size"
  jobMemory: "memory for this job, in Gb"
  minMemory: "Minimal amount of memory to assign to the task"
  timeout: "Timeout in hours, needed to override imposed limits"
}

Int allocatedMemory = if minMemory > round(jobMemory * scaleCoefficient) then minMemory else round(jobMemory * scaleCoefficient)

command <<<
 set -euxo pipefail
 ~{samtools} mpileup -q 1 -r ~{region} -f ~{refFasta} ~{inputNormal} ~{inputTumor} | awk -F "\t" '$4 > 0 && $7 > 0' | gzip -c > normtumor_sorted.pileup.gz 
>>>

runtime {
 memory: "~{allocatedMemory} GB"
 modules: "~{modules}"
 timeout: "~{timeout}"
}

output {
 File pileup = "normtumor_sorted.pileup.gz"
}
}

#=============================================================
# Task for concatenating CNV and SNV variants
#=============================================================
task mergeVariantsNative {
input {
 Array[File] filePaths
 String outputFile = "concatenated_variants"
 String outputExtension = "csv"
 Int jobMemory = 6
 Int timeout   = 10
}

parameter_meta {
  filePaths: "Array of pileup files to concatenate"
  jobMemory: "memory in GB for this job"
  outputExtension: "Extension of the output file"
  outputFile: "Name of the output file"
  timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
 set -euo pipefail
 head -n 1 ~{filePaths[0]} > "~{outputFile}.~{outputExtension}"
 cat ~{sep=' ' filePaths} | sort -V -k 1,2 | grep -v ^chrom | grep -v ^chrM >> "~{outputFile}.~{outputExtension}"
 cat ~{sep=' ' filePaths} | awk '{if($1 == "chrM"){print $0}}' | sort -V -k 1,2 >> "~{outputFile}.~{outputExtension}"
 if [ ! -s ~{outputFile}.~{outputExtension} ] ; then
  rm ~{outputFile}.~{outputExtension}
 fi
>>>

runtime {
 memory: "~{jobMemory} GB"
 timeout: "~{timeout}"
}

output {
  File mergedVariants = "~{outputFile}.~{outputExtension}"
}
}


#=============================================================
# Task for concatenating CNV and SNV variants
#=============================================================
task mergeVariantsVcf {
input {
 Array[File] filePaths
 String outputFile = "concatenated_vcf"
 String outputSuffix = "snp"
 String modules
 String seqDictionary 
 Int jobMemory = 12
 Int javaMemory = 8
 Int timeout   = 10
}

parameter_meta {
  filePaths: "Array of pileup files to concatenate"
  jobMemory: "memory in GB for this job"
  javaMemory: "memory in GB for java VM"
  outputFile: "Name of the output file"
  outputSuffix: "Suffix to use for an output file: snp or indel"
  seqDictionary: ".dict file for the reference in use"
  modules: "modules needed for this task"
  timeout: "Timeout in hours, needed to override imposed limits"
}

command<<<
 set -euxo pipefail
 unset _JAVA_OPTIONS
 java -Xmx~{javaMemory}G -jar $PICARD_ROOT/picard.jar SortVcf I=~{sep=' I=' filePaths} SD=~{seqDictionary} O=~{outputFile}.~{outputSuffix}.vcf
>>>

runtime {
 modules: "~{modules}"
 memory: "~{jobMemory} GB"
 timeout: "~{timeout}"
}

output {
  File mergedVcf = "~{outputFile}.~{outputSuffix}.vcf"
}
}

# ==========================================
#  configure and run Varscan in SNV mode
# ==========================================
task runVarscanSNV {
input {
  File inputPileup
  String sampleID ="VARSCAN"
  Float pValue = 0.05
  Int jobMemory  = 20
  Int minMemory = 4
  Int javaMemory = 6
  Int minCoverage = 8
  Int minCoverageNormal = 8
  Int minCoverageTumor = 6
  Float minVarFreq = 0.1
  Float minFreqForHom = 0.75
  Float normalPurity = 1.0
  Float tumorPurity = 1.0
  Float pValueHet = 0.99
  Float scaleCoefficient = 1.0
  Int strandFilter = 0
  Int validation = 0
  Int outputVcf = 0
  String logFile = "VARSCAN_SNV.log"
  String varScan = "$VARSCAN_ROOT/VarScan.jar"
  String modules = "varscan/2.4.2 java/8"
  Int timeout = 40
}

parameter_meta {
 inputPileup: "Input .pileup file for analysis"
 sampleID: "This is used as a prefix for output files"
 pValue: "somatic p-value for SNV calling, default is 0.05"
 minCoverage: "Minimum coverage in normal and tumor to call variant [8]"
 minCoverageNormal: "Minimum coverage in normal to call somatic [8]"
 minCoverageTumor: "Minimum coverage in tumor to call somatic [6]"
 minVarFreq: "Minimum variant frequency to call a heterozygote [0.10]"
 minFreqForHom: "Minimum frequency to call homozygote [0.75]"
 normalPurity: "Estimated purity (non-tumor content) of normal sample [1.00]"
 tumorPurity: "Estimated purity (tumor content) of normal sample [1.00]"
 pValueHet: "p-value threshold to call a heterozygote [0.99]"
 strandFilter: "If set to 1, removes variants with >90% strand bias"
 validation: "If set to 1, outputs all compared positions even if non-variant"
 scaleCoefficient: "Scaling coefficient for RAM allocation, depends on chromosome size"
 jobMemory: "Memory in Gb for this job"
 minMemory: "A minimum amount of memory allocated to the task, overrides the scaled RAM setting"
 javaMemory: "Memory in Gb for Java"
 logFile: "File for logging Varscan messages"
 outputVcf: "Flag that when set to 1 indicates that we need results in vcf format"
 varScan: "path to varscan .jar file"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
}

Int allocatedMemory = if minMemory > round(jobMemory * scaleCoefficient) then minMemory else round(jobMemory * scaleCoefficient)

command <<<
 unset _JAVA_OPTIONS
 set -euxo pipefail
 python<<CODE
 import os
 import re
 varscan = os.path.expandvars("~{varScan}")
 varscanCommand = "zcat ~{inputPileup} | java -Xmx~{javaMemory}G -jar " + varscan + " somatic -mpileup 1 --somatic-p-value ~{pValue}"

 if "~{minCoverageNormal}" != "8":
    varscanCommand += " --min-coverage-normal ~{minCoverageNormal}"
 if "~{minCoverageTumor}" != "6":
    varscanCommand += " --min-coverage-tumor ~{minCoverageTumor}"
 if "~{minVarFreq}" != "0.1":
    varscanCommand += " --min-var-freq ~{minVarFreq}"
 if "~{minFreqForHom}" != "0.75":
    varscanCommand += " --min-freq-for-hom ~{minFreqForHom}"
 if "~{normalPurity}" != "1.0":
    varscanCommand += " --normal-purity ~{normalPurity}"
 if "~{tumorPurity}" != "1.0":
    varscanCommand += " --tumor-purity ~{tumorPurity}"
 if "~{pValueHet}" != "0.99":
    varscanCommand += " --p-value ~{pValueHet}"
 if "~{strandFilter}" != "0":
    varscanCommand += " --strand-filter 1"
 if "~{validation}" != "0":
    varscanCommand += " --validation 1"
 if "~{outputVcf}" != "0":
    varscanCommand += " --output-vcf 1"
    varscanCommand += " --output-snp ~{sampleID}.snp.vcf --output-indel ~{sampleID}.indel.vcf"
 else:
    varscanCommand += " --output-snp ~{sampleID}.snp --output-indel ~{sampleID}.indel"

 cvg = ~{minCoverage}
 resultsOk = False
 f = open("~{logFile}", "w+")
 f.write('[Varscan log]' + '\n')

 while not resultsOk:
     run_command = varscanCommand
     run_command += " --min-coverage " + str(cvg)
     cvg -= 2
     f.write('[' + run_command + ']\n')
     res_string = os.popen(run_command + " 2>&1").read()
     f.write(res_string)
     m = re.search(r'(\d+)', str(res_string))
     m = re.search(r'(\d+) had sufficient coverage', res_string)
     if not m or (m and m.group(1) == '0' and cvg <= 2):
         f.write('Unable to run Varscan even with min-coverage set to ' + str(cvg) + ' aborting...\n')
         break
     elif m and m.group(1) != '0':
         resultsOk = True
         f.write('Success!\n')
     else:
         f.write('Coverage threshold too high, trying min-coverage ' + str(cvg) + '...\n')

 if not resultsOk:
     f.write('Varscan failed\n')
 f.close()
 CODE
>>>

runtime {
  memory: "~{allocatedMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File? snpFile   = "~{sampleID}.snp"
  File? indelFile = "~{sampleID}.indel"
  File? snpVcfFile = "~{sampleID}.snp.vcf"
  File? indelVcfFile = "~{sampleID}.indel.vcf"
}
}

# ==========================================
#  configure and run Varscan in CNV mode
# ==========================================
task runVarscanCNV {
input {
  File inputPileup
  String sampleID ="VARSCAN"
  Float pValue = 0.05
  Float scaleCoefficient = 1.0
  Int jobMemory  = 20
  Int minMemory = 4
  Int javaMemory = 6
  String logFile = "VARSCAN_CNV.log"
  String varScan = "$VARSCAN_ROOT/VarScan.jar"
  String modules = "varscan/2.4.2 java/8"
  Int timeout = 40
}

parameter_meta {
 inputPileup: "Input .pileup file for analysis"
 sampleID: "This is used as a prefix for output files"
 pValue: "p-value for cnv calling, default is 0.05"
 jobMemory: "Memory in Gb for this job"
 minMemory: "A minimum amount of memory allocated to the task, overrides the scaled RAM setting"
 javaMemory: "Memory in Gb for Java"
 logFile: "File for logging Varscan messages"
 scaleCoefficient: "Scaling coefficient for RAM allocation, depends on chromosome size"
 varScan: "path to varscan .jar file"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
}

Int allocatedMemory = if minMemory > round(jobMemory * scaleCoefficient) then minMemory else round(jobMemory * scaleCoefficient)

command <<<
 unset _JAVA_OPTIONS
 set -euxo pipefail
 python<<CODE
 import os
 import re
 varscan = os.path.expandvars("~{varScan}")
 varscanCommand = "zcat ~{inputPileup} | java -Xmx~{javaMemory}G -jar " + varscan + " copynumber --output-file ~{sampleID} -mpileup 1 --p-value ~{pValue}"
 cvg = 0
 resultsOk = False
 f = open("~{logFile}", "w+")
 f.write('[Varscan log]' + '\n')

 while not resultsOk:
     run_command = varscanCommand
     if cvg == 0:
         cvg = 19
     cvg -= 4
     run_command += " --min-coverage " + str(cvg)
     f.write('[' + run_command + ']\n')
     res_string = os.popen(run_command + " 2>&1").read()
     m = re.search(r'(\d+)', str(res_string))
     m = re.search(r'(\d+) had sufficient coverage', res_string)
     if not m or (m and m.group(1) == '0' and cvg <= 2):
         f.write('Unable to run Varscan even with min-coverage set to ' + str(cvg) + ' aborting...\n')
         break
     elif m and m.group(1) != '0':
         resultsOk = True
         f.write('Success!\n')
     else:
         f.write('Coverage threshold too high, trying min-coverage ' + str(cvg) + '...\n')

 if not resultsOk:
     f.write('Varscan failed\n')
 f.close()
 CODE
>>>

runtime {
  memory: "~{allocatedMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File? resultFile = "~{sampleID}.copynumber"
}
}

# ====================================================
#  Additional smoothing for Varscan data. Visualization
#  will be implemented after updates to rstats module
#       needed to enable png/jpg rendering
# ======================================================

task smoothData {
input {
 File copyNumberFile
 String varScan = "$VARSCAN_ROOT/VarScan.jar"
 String modules = "varscan/2.4.2 java/8 rstats/3.6"
 Int min_coverage  = 20
 Int max_homdel_coverage = 5
 Int min_tumor_coverage = 10
 Float del_threshold = 0.25
 Float amp_threshold = 0.25
 Int min_region_size = 10
 Int recenter_up = 0
 Int recenter_down = 0
 String sampleID ="VARSCAN"
 Int jobMemory  = 16
 Int javaMemory = 6
}

parameter_meta {
 copyNumberFile: "Output from Varscan"
 varScan: "Path to VarScan jar file"
 modules: "Modules for this job"
 min_coverage: "Fine-tuning parameter for VarScan"
 max_homdel_coverage: "Max coverage form homozygous deletion, default is 5"
 min_tumor_coverage: "Min coverage in tumor sample, default is 10"
 del_threshold: "Fine-tuning parameter for VarScan"
 amp_threshold: "Amplification threshold to report, default is 0.25"
 min_region_size: "Fine-tuning parameter for VarScan"
 recenter_up: "Fine-tuning parameter for VarScan"
 recenter_down: "Fine-tuning parameter for VarScan"
 sampleID: "sample id (used as prefix for result files)"
 jobMemory: "Memory in Gb for this job"
 javaMemory: "Memory in Gb for Java"
}

command <<<
 python<<CODE
 import os
 varscan = os.path.expandvars("~{varScan}")
 filterCommand = "java -Xmx~{javaMemory}G -jar " + varscan + " copyCaller ~{copyNumberFile} --output-file ~{sampleID}.copynumber.filtered"

 if "~{min_coverage}" != "20":
    filterCommand += " --min-coverage ~{min_coverage}"
 if "~{min_tumor_coverage}" != "10":
    filterCommand += " --min-tumor-coverage ~{min_tumor_coverage}"
 if "~{max_homdel_coverage}" != "5":
    filterCommand += "--max-homdel-coverage ~{max_homdel_coverage}"
 if "~{del_threshold}" != "0.25":
    filterCommand += " --del-threshold ~{del_threshold}"
 if "~{amp_threshold}" != "0.25":
    filterCommand += " --amp-threshold ~{amp_threshold}"
 if "~{min_region_size}" != "10":
    filterCommand += " --min-region-size ~{min_region_size}"
 if "~{recenter_up}" != "0":
    filterCommand += " --recenter-up ~{recenter_up}"
 if "~{recenter_down}" != "0":
    filterCommand += " --recenter-down ~{recenter_down}"

 message = os.popen(filterCommand + " 2>&1").read()
 CODE
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
 File? filteredData = "~{sampleID}.copynumber.filtered"
}

}

