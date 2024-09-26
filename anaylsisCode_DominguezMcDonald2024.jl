#===========================
Bradon R. McDonald
brmcdonald@wisc.edu
Kisat Lab
University of Wisconsin-Madison
2024
==========================#

inDir = "./brakenCountData_DominguezMcDonald2024"
classStats = "./totalReadcounts_DominguezMcDonald2024.tsv"
outDir = "./analysisOutput_DominguezMcDonald2024"
dataName = "DominguezMcDonald2024"
pathYaml = "./pathogenData.yaml"

include("./anaylsisCode_DominguezMcDonald2024_functions.jl")
include("./anaylsisCode_DominguezMcDonald2024_plots.jl")
pathogensDict = YAML.load(open(pathYaml))

###########################################################################
# Running analysis
###########################################################################

run(`mkdir -p $outDir`)

D01 = parseTaxonData(inDir)
annotateDepth!(D01,"./totalReadcounts_DominguezMcDonald2024.tsv")
annotateUnmapped!(D01,classStats)
annotateInput!(D01,"./cfDNA_concentration_DominguezMcDonald2024.tsv")
dfBKA_mp,dfBKB_mp,BK_A_mp,BK_B_mp = backgroundGenus(D01,pathogensDict,metric="microprop")

dfMaster = buildKnownTable(D01,["B251","B266","B268","B297","B304"],pathogensDict)
dfMasterPropMicro = buildKnownTable(D01,["B251","B266","B268","B297","B304"],pathogensDict,metric="propMicro")
dfMasterDepth = buildKnownTable(D01,["B251","B266","B268","B297","B304"],pathogensDict,metric="depth")
dfMasterCount = buildKnownTable(D01,["B251","B266","B268","B297","B304"],pathogensDict,metric="count")
dfMasterTotal = buildKnownTable(D01,["B251","B266","B268","B297","B304"],pathogensDict,metric="total")

dfMasterADM_mp = buildKnownTable(D01,["B251","B266","B268","B297","B304"],pathogensDict,metric="adm_mp",BK=BK_B_mp)

CSV.write("$(outDir)/positivePatientTable_$(dataName).tsv",dfMaster,delim='\t')
CSV.write("$(outDir)/positivePatientTable_$(dataName)_depth.tsv",dfMasterDepth,delim='\t')
CSV.write("$(outDir)/positivePatientTable_$(dataName)_counts.tsv",dfMasterCount,delim='\t')
CSV.write("$(outDir)/positivePatientTable_$(dataName)_totalMicrobial.tsv",dfMasterTotal,delim='\t')
CSV.write("$(outDir)/positivePatientTable_$(dataName)_ADM.tsv",dfMasterADM_mp,delim='\t')

dfTM_total = microFractionSummary(D01)
dfTM_unmapped = microFractionSummary(D01,denom="unmapped")

MI_01_SS = microFractionComparison(D01,prepType="SS")
MI_01_EMS = microFractionComparison(D01,prepType="EMS")
MI_01_STD = microFractionComparison(D01,prepType="STD")
MI_01_SSS = microFractionComparison(D01,prepType="SS_EMS")

#Boxplots of total microbial fraction
tmData_STD = vcat([MI_01_STD["mfA"],MI_01_STD["mfBPos"],MI_01_STD["mfBneg"]]...)
tmData_EMS = vcat([MI_01_EMS["mfA"],MI_01_EMS["mfBPos"],MI_01_EMS["mfBneg"]]...)
tmData_SS = vcat([MI_01_SS["mfA"],MI_01_SS["mfBPos"],MI_01_SS["mfBneg"]]...)
tmData_SSS = vcat([MI_01_SSS["mfA"],MI_01_SSS["mfBPos"],MI_01_SSS["mfBneg"]]...)

prepLabels = ["Takara","Takara\nSize selected","SRSLY","SRSLY\nSize selected"]
p0 = microFractionBoxplot2([tmData_STD,tmData_EMS,tmData_SS,tmData_SSS],prepLabels;title="Microbial Fragment Fraction");
save("$(outDir)/microbialFraction_$(dataName).pdf",p0)

#Comparing total microbial fraction vs cfDNA yield
YVT = Dict()
YVT["STD"] = [vcat([MI_01_STD["inputA"],MI_01_STD["inputBPos"],MI_01_STD["inputBneg"]]...),vcat([MI_01_STD["mfA"],MI_01_STD["mfBPos"],MI_01_STD["mfBneg"]]...)]
YVT["EMS"] = [vcat([MI_01_EMS["inputA"],MI_01_EMS["inputBPos"],MI_01_EMS["inputBneg"]]...),vcat([MI_01_EMS["mfA"],MI_01_EMS["mfBPos"],MI_01_EMS["mfBneg"]]...)]
YVT["SS"] = [vcat([MI_01_SS["inputA"],MI_01_SS["inputBPos"],MI_01_SS["inputBneg"]]...),vcat([MI_01_SS["mfA"],MI_01_SS["mfBPos"],MI_01_SS["mfBneg"]]...)]
YVT["SS_EMS"] = [vcat([MI_01_SSS["inputA"],MI_01_SSS["inputBPos"],MI_01_SSS["inputBneg"]]...),vcat([MI_01_SSS["mfA"],MI_01_SSS["mfBPos"],MI_01_SSS["mfBneg"]]...)]

p1 = scatterPlot(YVT["STD"][1],YVT["STD"][2],
	logY=true,logX=true,w=3.5,h=3.5,ptSize=8,tickSize=14,title="cfDNA vs Total Microbial Takara");
p2 = scatterPlot(YVT["EMS"][1],YVT["EMS"][2],
	logY=true,logX=true,w=3.5,h=3.5,ptSize=8,tickSize=14,title="cfDNA vs Total Microbial Takara-Size Selected");
p3 = scatterPlot(YVT["SS"][1],YVT["SS"][2],
	logY=true,logX=true,w=3.5,h=3.5,ptSize=8,tickSize=14,title="cfDNA vs Total Microbial SRSLY");
p4 = scatterPlot(YVT["SS_EMS"][1],YVT["SS_EMS"][2],
	logY=true,logX=true,w=3.5,h=3.5,ptSize=8,tickSize=14,title="cfDNA vs Total Microbial SRSLY-SS");

save("$(outDir)/cfDNA_vs_microbialFraction_std_$(dataName).pdf",p1)
save("$(outDir)/cfDNA_vs_microbialFraction_ems_$(dataName).pdf",p2)
save("$(outDir)/cfDNA_vs_microbialFraction_ss_$(dataName).pdf",p3)
save("$(outDir)/cfDNA_vs_microbialFraction_ss-ss_$(dataName).pdf",p4)

##############################
# Patient longitudinal plots
##############################
pats = ["B251","B266","B268","B297","B304"]
preps = ["STD","EMS","SS","SS_EMS"]
pathogens = ["Streptococcus","Staphylococcus","Klebsiella","Haemophilus","Staphylococcus"]
colorSet = ["blue","red","purple","teal"]

CTP = Dict{String,Vector{Int}}()
CTP["B251"] = [2,9]
CTP["B266"] = [1]
CTP["B268"] = [6]
CTP["B297"] = [7]
CTP["B304"] = [2]

CTN = Dict{String,Vector{Int}}()
CTN["B251"] = [1,3,8]
CTN["B266"] = [5]
CTN["B268"] = [1,5,10]
CTN["B297"] = [6,8]
CTN["B304"] = [1]


outF = "./$(outDir)/posPatients_Lon_PathYieldV4_sig_microProp"
plotPatientLonADM(dfMasterPropMicro,dfMasterADM_mp,["STD","EMS","SS","SS_EMS"],outF,yaxis="microprop")

###################################
#	Stats
###################################

dfM = CSV.File(classStats) |> DataFrame
dsDNA = filter(:SampleName=>x->occursin("SS-EMS",x) == false && occursin("TAK",x) == false && occursin("-SS",x) == false,dfM)
dsDNAss = filter(:SampleName=>x->occursin("SS-EMS",x) == false && occursin("-TAK",x),dfM)
ssDNA = filter(:SampleName=>x->occursin("SS-EMS",x) == false && occursin("-SS",x),dfM)
ssDNAss = filter(:SampleName=>x->occursin("SS-EMS",x),dfM)

#mean unmapped fraction in dsDNA vs ssDNA, enrichment fold
aData =[r.GoodUnmapped/r.GoodFragments for r in eachrow(dsDNA)]
bData = [r.GoodUnmapped/r.GoodFragments for r in eachrow(ssDNA)]
cData =[r.GoodUnmapped/r.GoodFragments for r in eachrow(dsDNAss)]
dData = [r.GoodUnmapped/r.GoodFragments for r in eachrow(ssDNAss)]
aMean = mean(aData)
bMean = mean(bData)
cMean = mean(cData)
dMean = mean(dData)
enrich = [(bData[i]-aData[i])/aData[i] for i in 1:length(aData)]

#mean classification prop in dsDNA and ssDNA, enrichment fold
aData = [sum(values(D01[x.SampleName].phylumData)) / D01[x.SampleName].totalFrags for x in eachrow(dsDNA)]
bData = [sum(values(D01[x.SampleName].phylumData)) / D01[x.SampleName].totalFrags for x in eachrow(ssDNA)]

mean(aData)
std(aData)
mean(bData)
std(bData)

enrich = [(bData[i]-aData[i])/aData[i] for i in 1:length(aData)]
mean(enrich)
std(enrich)

#Mean micro fragments in all preps, fold enrichment from dsDNA to ssDNA-SS
aData = [sum(values(D01[x.SampleName].phylumData)) / D01[x.SampleName].totalFrags for x in eachrow(dsDNA)]
bData = [sum(values(D01[x.SampleName].phylumData)) / D01[x.SampleName].totalFrags for x in eachrow(ssDNA)]
cData = [sum(values(D01[x.SampleName].phylumData)) / D01[x.SampleName].totalFrags for x in eachrow(dsDNAss)]
dData = [sum(values(D01[x.SampleName].phylumData)) / D01[x.SampleName].totalFrags for x in eachrow(ssDNAss)]

mean(aData)
mean(bData)
mean(cData)
mean(dData)

enrich = [(dData[i]-aData[i])/aData[i] for i in 1:length(aData)]

#Yield vs microbial fraction
allYield = [D01[x].cfdna for x in keys(D01) if D01[x].patient in Set(["B251","B266","B268","B297","B304"]) && D01[x].prepType == "STD"]
