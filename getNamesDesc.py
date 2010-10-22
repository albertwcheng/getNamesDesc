#!/usr/bin/python


'''



Copyright 2010 Wu Albert Cheng <albertwcheng@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

'''
"""

Get official gene names and descriptions, aliases, GO terms, etc.
Upated: 3/3/2009

updated 11/18 : --correct-excel

## Example 
## ./getNamesDesc.py --useProbeDictionary=/net/coldfact/data/awcheng/Jaenisch/CHRIS-MIcroarrays/Qiagen_32K_v3.0_master_new.probemap --useAlsoProbeCol1=1 mm9.config /net/coldfact/data/awcheng/Jaenisch/CHRIS-MIcroarrays/february-MSI-ARRAYS/Normalized\ data/NEWER_Analysis_all/MsiRG200.txt 2 3 > tryprobe.out 2> tryprobe.err2

"""
import sys;
import os.path;
from configobj import ConfigObj;
from getopt import getopt;
from albertcommon import *

watchList=[];




def getScriptLocation():
	fullpath=os.path.realpath(sys.argv[0]);
	(head,tail)=os.path.split(fullpath);
	return head;

def globalize(D):
	for k in D.keys():
		globals()[k]=D[k];


def linkPubMed(queryTerm):
	return "http://www.ncbi.nlm.nih.gov/sites/entrez?EntrezSystem2.PEntrez.Pubmed.Pubmed_SearchBar.Db=pubmed&EntrezSystem2.PEntrez.Pubmed.Pubmed_SearchBar.Term="+queryTerm#+"&EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_DisplayBar.sPresentation=Abstract&EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_DisplayBar.sSendTo=Text";


def getKeyWordWatchList(filename):
	global watchList;
	fil=open(filename);
	for line in fil:
		line=line.strip();
		if(line==""):
			continue;
		watchList.append(line.upper());
		

def parseEntryFromData(arr,datalines):
	curEntry=dict();
	for line in datalines:
		line=line.strip();
		if line=="":
			if len(curEntry.keys())>0	:
				arr.append(curEntry);
			curEntry=dict();
		else:
			key=line[0:2];	
			value=line[3:];
			if(key=="%R"):
				curEntry["PMID"]=value;
			elif(key=="%T"):
				curEntry["Title"]=value;
			elif(key=="%J"):
				curEntry["Journal"]=value;
			elif(key=="%V"):
				curEntry["Volume"]=value;
			elif(key=="%N"):
				curEntry["Number"]=value;
			elif(key=="%8"):
				curEntry["Date"]=value;
			elif(key=="%C"):
				curEntry["Institute"]=value;
			elif(key=="%X"):
				curEntry["Abstract"]=value;
			elif(key=="%A"):
				if(not curEntry.has_key("Authors")):
					curEntry["Authors"]=[];
				curEntry["Authors"].append(value);


def getGoPubMedEntry(geneName):
	arr=[];
	if not parseEntryFromFile(arr,	goPubMedPrefix+geneName+".endNote"):
		print >>sys.stderr,"No Annotation endNote for gene "+geneName;
	
	return arr;
		



def parseEntryFromFile(arr,filename):
	if not os.path.exists(filename):
		return False;

	fil=open(filename);
	
	parseEntryFromData(arr,fil);
	fil.close();
	return True;


def loadProbeDictionary(filename):
	print >>sys.stderr,"Getting Probe Dictionary from "+filename;
	fil=open(filename);
	ProbeDict=dict();

	for lin in fil:
		spliton=lin.strip().split("\t");
		if(len(spliton)<2):
			continue;
		probe=spliton[0].strip();
		geneID=spliton[1].strip();
		ProbeDict[probe]=geneID;
	fil.close();
	return ProbeDict;

def findFile(prefices,filename):
	#first without
	if os.path.exists(filename):
		return filename
	
	for prefix in prefices:
		newfilename=prefix+"/"+filename
		if os.path.exists(newfilename):
			return newfilename

	print >> sys.stderr,filename,"not found in",prefices
	raise IOError


def loadGO(GOT,filename):
	
	

	print >>sys.stderr,"Getting GO terms from "+filename;
	fil=open(filename);
	for lin in fil:
		spliton=lin.strip().split("\t");
		if(len(spliton)<2):
			continue;
		GOID=spliton[0].strip();
		GOTerm=spliton[1].strip();
		GOT[GOID]=GOTerm;
	fil.close();
	

def loadGOA(nameMaps,GOs,GOT,filename,GOName):
	print >> sys.stderr,"Getting GO annotation " + GOName+ " from "+filename;
	fil=open(filename);
	for lin in fil:
		spliton=lin.strip().split("\t");
		if(len(spliton)<3):
			continue;
		ensGeneID=spliton[0] #or a symbol .upper();
		GOEvi=spliton[2];
		GOID=spliton[1];
		if ensGeneID not in nameMaps:
			#print >> sys.stderr, ensGeneID +" not found, ignored";
			continue;
		
		#print >> sys.stderr, ensGeneID + "found";
		key=nameMaps[ensGeneID];
		if(len(GOID.strip())>0):
			if(GOID in GOT):
				GOID+=":"+GOT[GOID];
			GOs[key][GOName].add(GOID);
	print >>sys.stderr, "Done loading GO set "+GOName+" from "+filename;


def getGeneID(nameMaps,K):
	#K=K.upper();
	if nameMaps.has_key(K):
		return nameMaps[K];
	else:
		return "";


def correctExcel(searchTerms):
	additor=""
	for term in searchTerms:
		splits=term.strip().split("-")
		if len(splits)==2:
			if splits[1].upper()=="SEP":
				additor= "Sept"+splits[0]
				if not additor in searchTerms and not additor.upper() in searchTerms:
					searchTerms.append(additor)
					searchTerms.append(additor.upper())
			elif splits[1].upper()=="MAR":
				additor=  "March"+splits[0]
				if not additor in searchTerms and not additor.upper() in searchTerms:
					searchTerms.append(additor)
					searchTerms.append(additor.upper())				

			



def loadSymbolNamesDescExtra(nameMaps,geneInfo,filename):
	print >> sys.stderr,"Loading extra gene infos from "+filename;
	fil=open(filename);
	for lin in fil:

		#lin=lin.strip();
		spliton=lin.split("\t");
		if(len(spliton)<3):
			print >> sys.stderr, "line columns < 3, ignored:" + lin;
			continue;
		#print >> sys.stderr, spliton
		geneID=spliton[0];
		geneAliases=spliton[2].strip();
		geneAliasesSplit=geneAliases.split("|");
		geneAccessions=spliton[1];
		geneAccessionsSplit=geneAccessions.split("|");
		



		if(geneID==""):
			for g in geneAliasesSplit:
				geneID=getGeneID(nameMaps,g);


				if len(geneID)>0:

					break;
			
			if len(geneID)<1:
				geneID=geneAccessionsSplit[0];
				if geneID=="":
					continue;
				#print >> sys.stderr, "Gene Not Found, use accession=",geneID;
			#else:
			#	print >> sys.stderr, "Gene ID Found for",geneAccessionsSplit[0],"=",geneID;
		
		if not geneInfo.has_key(geneID):			
			
			geneInfo[geneID]=dict();
			geneInfo[geneID]["Desc"]="";
			geneInfo[geneID]["geneID"]=geneID;
			geneInfo[geneID]["Sym"]=geneAliasesSplit[0];
			geneInfo[geneID]["Aliases"]="";
			geneInfo[geneID]["Names"]=set();
			geneInfo[geneID]["Accessions"]="";
		

		#geneInfo[geneID]["Sym"]+=geneAliasesSplit[0];
		geneInfo[geneID]["Aliases"]+="|"+geneAliases;	

		for g in geneAliasesSplit:
			if g.strip()=="":
				continue;
			geneInfo[geneID]["Names"].add(g);

			nameMaps[g]=geneID; #Restored! .upper()
		geneInfo[geneID]["Accessions"]+="|"+geneAccessions;
		
		for g in geneAccessionsSplit:
			if g.strip()=="":
				continue
			gs=g.split(":");
			if(len(gs)<2):

				nameMaps[g]=geneID; #.upper()

				if (not geneInfo[geneID].has_key("Ensembl")) and g.find("ENS")>-1:
					geneInfo[geneID]["Ensembl"]=g;
					nameMaps[g]=geneID; #.upper()

				continue; ##???

			accType=gs[0];
			accID=gs[1].strip();
			if( len(accType)==0 or len(accID)==0 ):
				continue;
			geneInfo[geneID][accType]=accID;
			if(accType=="Ensembl"):

				nameMaps[accID]=geneID; #.upper()
		geneInfo[geneID]["GO"]=set();


		
		
	print >> sys.stderr, "finish loading NCBI gene_info";
	fil.close();
	

def loadSymbolNamesDescNCBI(nameMaps,geneInfo,filename):
	print >> sys.stderr,"Loading NCBI gene_info from " + filename;
	fil=open(filename);
	for lin in fil:
		lin=lin.strip();
		spliton=lin.split("\t");
		if(len(spliton)<14):
			print >> sys.stderr, "line columns < 14, ignored:" + lin;
			continue;
		
		geneID=spliton[1];
		geneSym=spliton[2];
		geneLocusTag=spliton[3]; #added	4/14/2009	
		geneAliases=spliton[4];
		geneAccessions=spliton[5];
		geneDesc=spliton[8];
		geneInfo[geneID]=dict();
		geneInfo[geneID]["Desc"]=geneDesc;
		geneInfo[geneID]["geneID"]=geneID;
		geneInfo[geneID]["Sym"]=geneSym;
		
		geneInfo[geneID]["Aliases"]=geneAliases;
		geneAliasesSplit=geneAliases.split("|");  
		geneAliasesSplit.append(geneLocusTag)#added	4/14/2009	

		geneInfo[geneID]["Names"]=set();
		geneInfo[geneID]["Names"].add(geneSym);
		nameMaps[geneSym]=geneID;# .upper()
		for g in geneAliasesSplit:
			if g.strip()=="":
				continue;
			geneInfo[geneID]["Names"].add(g);
			nameMaps[g]=geneID; #Restored!# .upper()
		geneInfo[geneID]["Accessions"]=geneAccessions;
		geneAccessionsSplit=geneAccessions.split("|");
		for g in geneAccessionsSplit:
			gs=g.split(":");
			if(len(gs)<2):
				continue;
			accType=gs[0];
			accID=gs[1];
			if( len(accType)==0 or len(accID)==0 ):
				continue;
			geneInfo[geneID][accType]=accID;
			if(accType=="Ensembl"):
				nameMaps[accID]=geneID; #.upper()
		geneInfo[geneID]["GO"]=set();
				
		
	print >> sys.stderr, "finish loading NCBI gene_info";
	fil.close();

def loadSymbolNamesDescEns(nameMaps,GOs,filename):
	print >> sys.stderr,"Loading reference ...";
	fil=open(filename);
	for lin in fil:
		lin=lin.strip();
		if len(lin)<1:
			continue;
		spliton=lin.split("\t");
		if len(spliton)<1:		
			continue;		
		ensGeneID=spliton[0];

		if(len(spliton)>=2):
			assGeneName=spliton[1];
		else:	
			assGeneName="";

		if(len(spliton)>=3):
			desc=spliton[2];
		else:
			desc="";
		if(len(spliton)>=5):
			UCSCID=spliton[4];
		else:	
			UCSCID="";
		if(len(spliton)>=4):
			entrezGeneID=spliton[3];
		else:
			entrezGeneID="";

		key=ensGeneID; ##+"/"+assGeneName;
		if not GOs.has_key(key):
			
			GOs[key]=dict();
			GOs[key]["ensGeneID"]=ensGeneID;
			GOs[key]["assGeneName"]=assGeneName;
			GOs[key]["UCSCID"]=set();
			GOs[key]["entrezGeneID"]=entrezGeneID;
			GOs[key]["aliases"]=set();
			GOs[key]["GO"]=set();
			GOs[key]["desc"]=desc;		
			
			#print >> sys.stderr, "add: ",
			#print >> sys.stderr,GOs[key];			
						
			nameMaps[ensGeneID]=key; #.upper()
			nameMaps[assGeneName]=key; #.upper()
			nameMaps[UCSCID]=key; #.upper()
			nameMaps[entrezGeneID]=key; #.upper()
		GOs[key]["UCSCID"].add(UCSCID);
		
							
	fil.close();
	print >> sys.stderr,"Done Loading Reference";




def main(filename,colGene1,startRow1,keywordWatchList,namesout,goPubMed,ProbeDict,useAlsoProbeCol1,onlyReplace,correctExcelFlag,replaceFirstOnly,printAllValuesOf,prefixLoading):
	namesout=open(namesout,"w");

	scriptPath=getScriptLocation()
	prefixLoading.append(scriptPath) 

	global taxID;
	global geneNameDescFileNCBI;
	global geneBGOFile;
	global geneGOFullFile;
	global watchList;
	global geneNameDescFileExtra;
	colGene=colGene1-1
	useAlsoProbeCol=useAlsoProbeCol1-1;
	GOT=dict();
	GOs=dict();
	nameMaps=dict();

	#load config
	geneNameDescFileNCBI=findFile(prefixLoading,geneNameDescFileNCBI)
	loadSymbolNamesDescNCBI(nameMaps,GOs,geneNameDescFileNCBI);

	for geneExtra in geneNameDescFileExtra:
		geneExtra=findFile(prefixLoading,geneExtra)
		loadSymbolNamesDescExtra(nameMaps,GOs,geneExtra)

	for geneGOFullF in geneGOFullFile:
		geneGOFullF=findFile(prefixLoading,geneGOFullF)
		loadGO(GOT,geneGOFullF);

	for geneBGOAF in geneBGOAFile:
		geneBGOAF=findFile(prefixLoading,geneBGOAF)
		loadGOA(nameMaps,GOs,GOT,geneBGOAF,"GO");

	getKeyWordWatchList(keywordWatchList);



	if printAllValuesOf!="":
		for GOKey,GOItem in GOs.items():
			printV=[GOKey]
			for v in printAllValuesOf:
				printV.append(GOItem[v])

			print >> sys.stdout,"\t".join(printV)
		return

	lino=0;	
	fil=generic_istream(filename);
	for lin in fil:
		lino+=1;
		lin=lin.rstrip("\r\n");
		if lino<startRow1:
			
			if not onlyReplace: #not just replace

				fieldsOut=[lin, "Defined Loci Gene #","Defined Loci Name","ensGeneID","geneSym","aliases","desc"]
				

				fieldsOut.append("iHOP")
				for term in watchList:
					fieldsOut.append( "GO:"+term)

				for term in watchList:
					fieldsOut.append(  "LIT:"+term)
				for term in watchList:
					fieldsOut.append(  term)
				fieldsOut.append(  "GO Full")

				if goPubMed:
					fieldsOut.append(  "GoPubMed...")
	
				print >> sys.stderr,"line 1 ignored";
				print >> sys.stdout, "\t".join(fieldsOut)

			else:
				print >> sys.stdout, lin

			continue
		
		spliton=lin.split("\t");
		
		if len(spliton) < colGene+1 and len(spliton)<useAlsoProbeCol+1:
			print >> sys.stderr, "spliton<1, ignored: line=",lino;
			print >> sys.stdout, lin
			continue;
		searchTerms=[];

		if len(spliton) >= colGene+1:
			geneIDString=spliton[colGene]
			if "///" in geneIDString:			
				searchTerms=geneIDString.split("///");
			else:
				searchTerms=geneIDString.split(",");
			

		if useAlsoProbeCol>=0:
			searchTerms.append(spliton[useAlsoProbeCol]);
		
		#print >>sys.stderr, "***********searchTerms: ",
		#print >>sys.stderr, searchTerms;
		
		###added nov 2009
		if correctExcelFlag:
			correctExcel(searchTerms)
		####

		keysForOutput=set();	
		
		for searchTerm in searchTerms:
			searchTerm=searchTerm.strip()
			if searchTerm=='-':
				#clean this stick case
				continue

			#print >> sys.stderr, "searchTerm loop:"+searchTerm;
			if ProbeDict.has_key(searchTerm):
				searchTerm=ProbeDict[searchTerm];
			
			searchTermUpper=searchTerm.strip(); #.upper()
			if len(searchTermUpper)==0:
				continue;
			if nameMaps.has_key(searchTermUpper):		
				keysForOutput.add(nameMaps[searchTermUpper]);

			comNameSplit=searchTermUpper.split("-"); #### removed in v2
			if(len(comNameSplit)>1):
				if comNameSplit[1]=="PENDING":				
					if nameMaps.has_key(comNameSplit[0]):
					#print >> sys.stderr, "adding "+comNameSplit[0];
						keysForOutput.add(nameMaps[comNameSplit[0]]);
		if not onlyReplace:
			fieldsOut=[ lin ]
		

		if(len(keysForOutput)==0):
			if onlyReplace:
				print >> sys.stdout, lin;		
			else:
				printBuf=str(len(keysForOutput))

			if len(spliton)>colGene:
				namesout.write(spliton[colGene]+"\n"); #still you need to add to the namesout, update other file as well
				print >> sys.stderr, str(len(keysForOutput))+" record(s) for "+spliton[colGene];
				if not onlyReplace: #added 4/14/2009
					printBuf+="\t"+spliton[colGene];
			else: #even out of range
				print >> sys.stderr, str(len(keysForOutput))+" record(s) for line "+lin;
			
			if not onlyReplace:
				fieldsOut.append( printBuf)
				print >> sys.stdout,"\t".join(fieldsOut)
				continue

		if(len(keysForOutput)>0):
			#print >> sys.stdout, "\t",
			#print >> sys.stdout, "\t",
			#print >> sys.stdout, "\t",
			#print >> sys.stdout,  "\t",
			#print >> sys.stdout, "\t",
			#print >> sys.stdout, "\t",
			#print >> sys.stdout,  "\t",
			#print >> sys.stdout, "\t",
		
		#else:	

			if(len(keysForOutput)>1):
				newkeysForOutput=[];
				for key in keysForOutput:
					matched=False;					
					GOItem=GOs[key];
					for searchTerm in searchTerms:
						if searchTerm in GOItem["Names"]:
							matched=True;

					if matched:
						newkeysForOutput.append(key);

				if len(newkeysForOutput)==1:
					keysForOutput=newkeysForOutput;

				elif len(newkeysForOutput)>1:
					newsymKeys=[];
					for key in keysForOutput:
						matched=False;					
						GOItem=GOs[key];
						for searchTerm in searchTerms:
							if searchTerm==GOItem["Sym"]:
								matched=True;

						if matched:
							newsymKeys.append(key);
					
					if len(newsymKeys)==1:
						keysForOutput=newsymKeys;
					

	

			arrSym=[];
			for key in keysForOutput:
				arrSym.append(GOs[key]["Sym"]);


			
			if onlyReplace:
				symsReplace=[];
				for key in keysForOutput:
					try:
						symsReplace.append(GOs[key][onlyReplace])
					except KeyError:
						print >> sys.stderr,"item",GOs[key]["Sym"],"does not have info for",onlyReplace
						#print >> sys.stderr, GOs
						

				if(len(symsReplace)>0):
					replaceString="|".join(symsReplace)
					if replaceFirstOnly:
						replaceString=symsReplace[0]
					if replaceString!="": #added 4/14/2009
						spliton[colGene]=replaceString

				#fieldsOut=spliton[:]  #fieldsOut.append( "\t".join(spliton))
				print >> sys.stdout,"\t".join(spliton)
			else:
				fieldsOut.append( str(len(keysForOutput)) )
				
				commSym=",".join(arrSym)
				if commSym=="": #added 4/14/2009
					commSym=spliton[colGene]

				fieldsOut.append( commSym )

				for key in keysForOutput:
					GOItem=GOs[key];
					if(len(GOItem["Sym"])>0):
						namesout.write(GOItem["Sym"]+"\n");
					if("Ensembl" in GOItem):
						fieldsOut.append(  GOItem["Ensembl"])
					else:
						fieldsOut.append( "")

					fieldsOut.append( GOItem["Sym"] )
					fieldsOut.append(  "|".join(GOItem["Names"]))
					fieldsOut.append(  GOItem["Desc"])
					fieldsOut.append(  "=HYPERLINK(\"http://www.ihop-net.org/UniPub/iHOP/in?ncbi_tax_id_1="+taxID+"&syns_1="+GOItem["Sym"]+"\",\"iHOP\")") ####
				
					strGOsSearchee="";
					strGOs="";
					strGOs+= "|".join(GOItem["GO"])+"\t";
					strGOsSearchee=strGOs;
					strLitSearchee="";
					
					if goPubMed:
						arrGPM=getGoPubMedEntry(GOItem["Sym"]);
						for lit in arrGPM:
							if goPubMed:
								strGOs+="=HYPERLINK(\""+linkPubMed(lit["PMID"])+"\",\""+lit["Title"]+"\")\t"+lit["Abstract"]+"\t";
					
							strLitSearchee+=lit["Title"]+" " +lit["Abstract"]+" ";
				
					strGOsSearcheeUpper=strGOsSearchee.upper();
					strLitSearcheeUpper=strLitSearchee.upper();
				

					sumCount=[];
					for term in watchList:
						#count term in strGOs
						count=strGOsSearcheeUpper.count(term);
						sumCount.append(count);		
						fieldsOut.append( str(count) )
					i=0
					for term in watchList:
						#count term in strGOs
						count=strLitSearcheeUpper.count(term);
					
						sumCount[i]+=count;		
						i+=1;					
						fieldsOut.append(str(count))

					for count in sumCount:
						#total
						fieldsOut.append(str(count))	
			
					fieldsOut.append( strGOs.strip())

				print >> sys.stdout,"\t".join(fieldsOut)
		
	print >> sys.stderr,"<Done>";
	fil.close();
	namesout.close();




if len(sys.argv)<5:
	print >>sys.stderr,"Usage: "+sys.argv[0]+" configName filename colGene1 startRow1"; ##########
	print >>sys.stderr,"Options:";
	print >>sys.stderr,"--nameout=";
	print >>sys.stderr,"--keyword-watch-list=";
	print >>sys.stderr,"--goPubMed"
	print >>sys.stderr,"--onlyReplace"
	print >>sys.stderr,"--headerRow"
	print >> sys.stderr,"--correct-excel"
	print >> sys.stderr,"--print-all-values-of"
	explainColumns(sys.stderr)
	exit();

opts,args=getopt(sys.argv[1:],'',['nameout=','keyword-watch-list=','goPubMed','useProbeDictionary=','useAlsoProbeCol1=','onlyReplace=','onlyReplaceFirstOf=','headerRow=','correct-excel','print-all-values-of=']);
print >> sys.stderr, args;

nameout="/dev/null";
watchlist="/dev/null";
goPubMed=False;
onlyReplace="";
fs="\t"
ProbeDict=dict();
useAlsoProbeCol1=0;
headerRow=-1
replaceFirstOnly=False
correctExcelFlag=False
printAllValuesOf=""

for o,a in opts:
	if o=='--nameout':
		nameout=a;
	elif o=='--keyword-watch-list':
		watchlist=a;
	elif o=='--goPubMed':
		goPubMed=True;
	elif o=='--useProbeDictionary':
		ProbeDict=loadProbeDictionary(a);
	elif o=='--useAlsoProbeCol1':
		useAlsoProbeCol1=int(a);
	elif o=='--onlyReplace':
		onlyReplace=a
	elif o=='--onlyReplaceFirstOf':
		replaceFirstOnly=True
		onlyReplace=a
	elif o=='--headerRow':
		headerRow=int(a)
	elif o=='--correct-excel':
		correctExcelFlag=True
	elif o=="--print-all-values-of":
		printAllValuesOf=a.split(",")



configName,filename,colGene1,startRow1=args;

startRow1=int(startRow1)

if headerRow==-1:
	headerRow=startRow1-1
	

try:
	colGene1=int(colGene1)
except ValueError:
	#not just a number: need something special
	header,prestarts=getHeader(filename,headerRow,startRow1,fs)
	colGene1=getCol0ListFromCol1ListStringAdv(header,colGene1)[0]+1
	print >> sys.stderr,"colGene1 is",colGene1

if not os.path.isfile(configName):
	scriptPath=getScriptLocation();
	print >> sys.stderr, configName,"not found in working directory. Using scriptPath",scriptPath;
	newName=scriptPath+"/"+configName#"mm9.config";	
	if os.path.isfile(newName):
		configName=newName;
	else:
		print >> sys.stderr, "Config file not found, exiting";
		exit();
	

config=ConfigObj(configName);

if "prefix" not in config:
	prefix=[""]

globalize(config);



#filename,colGene1,startRow1,keywordWatchList,namesout

main(filename,colGene1,startRow1,watchlist,nameout,goPubMed,ProbeDict,useAlsoProbeCol1,onlyReplace,correctExcelFlag,replaceFirstOnly,printAllValuesOf,prefix);

