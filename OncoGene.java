package edu.archive;

import org.json.simple.parser.ParseException;
import org.xml.sax.SAXException;

import javax.xml.parsers.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class OncoGene {
	
	public static ArrayList<String> in = new ArrayList<>();
	public static ArrayList<String> out = new ArrayList<>();
	public static HashMap<String,String> id = new HashMap<String,String>();
	public static HashMap<String,HashMap<String,String>> FH = new HashMap<String,HashMap<String,String>>();
	public static String[] stopwords = {"and", "i", "me", "my", "myself", "we", "our", "ours", "ourselves", "you", "your", "yours", "yourself", "yourselves", "he", "him", "his", "himself", "she", "her", "hers", "herself", "it", "its", "itself", "they", "them", "their", "theirs", "themselves", "what", "which", "who", "whom", "this", "that", "these", "those", "am", "is", "are", "was", "were", "be", "been", "being", "have", "has", "had", "having", "do", "does", "did", "doing", "a", "an", "the", "but", "if", "because", "as", "until", "while", "of", "at", "by", "for", "with", "about", "against", "between", "into", "through", "during", "before", "after", "above", "below", "to", "from", "up", "down", "in", "out", "on", "off", "over", "under", "again", "further", "then", "once", "here", "there", "when", "where", "why", "how", "all", "any", "both", "each", "few", "more", "most", "other", "some", "such", "nor", "only", "own", "same", "so", "than", "too", "very", "can", "will", "just", "don", "should", "now"};
	public static String folder = "/Users/m179121/Downloads/";
	
    public static void main(String[] args) throws IOException, ParseException, ParserConfigurationException, SAXException, java.text.ParseException {
    	
    	getOncoGene(folder+"in.txt");
    }
    
    public static void getOncoGene(String file) throws FileNotFoundException, IOException {
		
	    try (BufferedReader br = new BufferedReader(new FileReader(file))) {
	
	        String line = ""; 
	
	        while ((line = br.readLine()) != null) {
	        	
	        	String topic = ".";
	        	line = line.replaceAll("BRCA1 and 2|BRCA1 or BRCA2|BRCA1 and BRCA2|BRCA1/BRCA2|BRCA1 or 2", "BRCA1/2");
	        	
	        	if(line.startsWith("PATIENT")) {
	        		out.add("GENE\tANNOTATION\tTOPIC\tMOD\tVAR");
	        		continue;
	        	}
	        	
//	        	String gene = line.split("\t")[7];
//	        	if(gene.equals("met")||gene.equals("Met")) {
//	        		out.add("");
//	        		continue;
//	        	}
//	        	String sentence = line.split("\t")[8];
//	        	for (int i = 9; i < line.split("\t").length; i++) {
//	        		sentence = sentence + " " + line.split("\t")[i];
//	        	}
//
//	        	if(!sentence.replace(gene+" inhibitor", "inhibitor").contains(gene)) {
//	        		out.add("");
//	        		continue;
//	        	}
	        	
//	        	String gene = "BRCA1";
//	        	String sentence = "Of note, she has a history of breast cancer, and since she was last seen, we have received results of her mutation testing, and there are no mutations identified in 25 genes associated with hereditary cancer predispositions which included BARD1, BRCA1, BRCA2, CHEK2, RAD51C, RAD51D and TP53.";
	        	
	        	System.out.println(line);
	        	String gene = line.split("\t")[0];
	        	String sentence = line.split("\t")[1];
	        	
	        	sentence = sentence.replaceAll("(?i)(variant )?(of )?un[a-z\\s]+significance", "vus").replaceAll("(?i)( neg |not have|not has|wild type|wild-type|wild)", " negative ");
	        	
	        	if(gene.toLowerCase().equals("tp53")) {
	        		sentence = sentence.replaceAll("(?i)p53", gene).replaceAll("(?i)t"+gene, gene);
        		}
	        	if(gene.toLowerCase().equals("erbb2")) {
	        		sentence = sentence.replaceAll("(?i)her2", gene);
        		}
	        	if(gene.toLowerCase().equals("her2")) {
	        		sentence = sentence.replaceAll("(?i)erbb2", gene);
        		}
	        	
	        	sentence = removeStopWords(sentence).replaceAll("[\\p{Punct}&&[^$\\|\\./*>]]+", " ");
	        	if(sentence.endsWith(".")) {
	        		sentence = sentence.substring(0, sentence.length()-1);
	        	}

	        	Pattern pat = Pattern.compile("\\b(?i)(([A-Z][0-9]{2,}[A-Z](>[A-Z])?)|c\\.[a-z0-9>_]{1,}|p\\.[a-z0-9*>_]{1,}|loss|frameshift|translocat|splice( site)?( [ATGC>0-9]{1,})?|rearrange|insert|delet|amplification|duplication|fusion)[a-z]{0,}\\b");
	        	Matcher mat = pat.matcher(sentence);
	        	Pattern pat1 = Pattern.compile("\\b(?i)(mutation|alteration|variant)(s)?\\b");
	        	Matcher mat1 = pat1.matcher(sentence);
	        	Pattern pat2 = Pattern.compile("\\b(RAD51C|RAD51D|ABL1|ACVR1B|AKT1|AKT2|AKT3|ALK|APC|AR|ARAF|ARFRP1|ARID1A|ARID1B|ARID2|ASXL1|ATM|ATR|ATRX|AURKA|AURKB|AXIN1|AXL|B2M|BAP1|BARD1|BCL2|BCL2L1|BCL2L2|BCL6|BCL7A|BCOR|BCORL1|BIRC3|BLM|BRAF|BRCA1/2|BRCA1|BRCA2|BRD4|BRIP1|BRSK1|BTG1|BTG2|C11orf30|C17orf39|CALR|CAMTA1|CARD11|CBL|CCND1|CCND2|CCND3|CCNE1|CCT6B|CD274|CD36|CD58|CD70|CDC73|CDH1|CDK12| CDK4/6|CDK4|CDK6|CDK8|CDKN1A|CDKN1B|CDKN2A/B|CDKN2A|CDKN2B|CDKN2C|CEBPA|CHD2|CHD4|CHEK2|CIC|CIITA|CKS1B|CPS1|CREBBP|CRKL|CSF1R|CTCF|CTNNA1|CTNNB1|CUL3|CUX1|CXCR4|CYLD|DAXX|DDR2|DDX3X|DICER1|DNMT3A|EED|EGFR|EMSY|EP300|EPHA3|EPHA5|EPHA7|EPHB1|HER2|ERBB2|ERBB3|ERBB4|ERRFI1|ESR1|ETV1|ETV6|EWSR1|EZH2|FAF1|FAM123B|FANCA|FANCC|FANCD2|FANCE|FANCF|FANCG|FANCL|FAS|FAT1|FBXW7|FGF10|FGF14|FGF19|FGF23|FGF3|FGF4|FGF6|FGFR1|FGFR2|FGFR3|FGFR4|FH|FLCN|FLT1|FLT3|FLT4|FOXL2|FOXP1|FRS2|FUBP1|GABRA6|GATA1|GATA2|GATA3|GATA4|GATA6|GLI1|GNA11|GNA12|GNAQ|GNAS|GRIN2A|GRM3|GSK3B|H3F3A|HEY1|HGF|HIST1H1C|HIST1H1D|HIST1H2AM|HMGA2|HNF1A|HRAS|ID3|IDH1|IDH2|IGF1R|IGF2|IGH|IGK|IKBKE|IKZF3|IL7R|INHBA|INPP4B|IRF2|IRF4|IRS2|JAK1|JAK2|JAK3|JAZF1|JUN|KDM4C|KDM5A|KDM5C|KDM6A|KDR|KEAP1|KEL|KIT|KLHL6|KMT2C |KRAS|LEF1|LPP|LRP1B|LYN|LZTR1|MAGI2|MALT1|MAP2K1|MAP2K1 |MAP2K2 |MAP2K4|MAP3K1|MAP3K14|MAP3K6|MCL1|MDM2|MDM4|MED12|MEF2B|MEN1|MET|MITF|MLH1|MLL|MLL2|MLL3|MLLT10|MPL|MRE11A|MSH2|MSH6|MTOR|MUTYH|MYB|MYC|MYCL1|MYCN|MYD88|MYST3|NCOR2|NF1|NF2|NFE2L2|NFKBIA|NKX2-1|NOD1|NOTCH1|NOTCH2|NOTCH3|NPM1|NRAS|NSD1|NTRK1|NTRK3|NUP93|NUP98|PALB2|PARK2|PASK|PAX3|PAX5|PBRM1|PDCD1LG2 |PDGFRA|PDGFRB|PHF6|PICALM|PIK3C2B|PIK3CA|PIK3CB|PIK3CG|PIK3R1|PIK3R2|PIM1|PLCG2|PMS2|POLD1|POLE|POT1|PPP2R1A|PRDM1|PREX2|PRKAR1A|PRKCI|PRKDC|PRSS8|PTCH1|PTEN|PTPN11|PTPN2|PTPN6|QKI|RAC1|RAD21|RAD50|RAD51|RAF1|RANBP2|RARA|RB1|RBM10|RET|RICTOR|RNF43|ROS1|RPTOR|RUNX1|RUNX1T1|SDHA|SDHB|SETD2|SF3B1|SGK1|SLIT2|SMAD2|SMAD3|SMAD4|SMARCA4|SMARCB1|SMO|SNCAIP|SOCS1|SOCS3|SOX2|SOX9|SPEN|SPOP|SPTA1|SRC|SRSF2|SS18|STAG2|STAT3|STAT4|STAT6|STK11|SUFU|SUZ12|TAF1|TAF15|TBL1XR1|TBX3|TERC|TERT|TET2|TGFBR2|TMEM30A|TMPRSS2|TNFAIP3|TNFRSF11A|TNFRSF14|TOP1|TOP2A|TP53|TP63|TRAF2|TSC1|TSC2|TSHR|TYK2|U2AF1|VEGFA|VHL|WHSC1|WT1|XPO1|ZBTB2|ZMYM3|ZNF217|ZNF703)\\b");
	        	Matcher mat2 = pat2.matcher(sentence);
	        	Pattern pat3 = Pattern.compile("\\(.*\\)");
	        	Matcher mat3 = pat3.matcher(sentence);
	        	Pattern pat4 = Pattern.compile("\\b(?i)(pursue|elect|send|sent|pending|order|submit|submitted)(ing|es|ed|le|ility|e|s|d)?\\b");
	        	Pattern pat5 = Pattern.compile("\\b(?i)(analyze|review|discuss|recommend|encourag|explain|women|individual|male|female|lifetime)(ing|es|ed|le|ility|e|s|d)?\\b");
	        	Pattern pat6 = Pattern.compile("\\b(?i)(panel|Jewish|HBOC|autosomal|etiology|pathway|NCCN|includ)(ing|es|ed|le|ility|e|s|d)?\\b"); 
	        	Pattern pat7 = Pattern.compile("\\b(?i)(underlying|evaluation|personal|table|model|odds|probability|estimat)(ing|es|ed|le|ility|e|s|d)?\\b");
	        	Pattern pat8 = Pattern.compile("\\b(?i)(father|mother|parent|son|sister|sibling|daughter|child|cousin|brother|niece|nephew|uncle|aunt|grandparent|grandmother|grandfather|paternal|maternal|relative|family [a-z\\s]{0,20}history)(s)?\\b");
	        	Pattern pat9 = Pattern.compile("\\b(?i)(negative|(not|no|none) [a-z\\s]{0,15}pathogenic|not|no|none|normal|known|pathogenic|positive|vus|un[a-z]{0,10} significance|might|would|if)\\b");
	        	Pattern pat10 = Pattern.compile("\\b(?i)(approximately|hereditary|risk|suggestive|susceptib|likelihood|screen|chance|possib|associated)(ing|es|ed|le|ility|e|s|d)?\\b");
	        			
	        	while(mat2.find()) {//get gene name
//	        		System.out.println(mat2.group());
	        		if(!mat2.group().equals(gene)) {
	        			sentence = sentence.replace(mat2.group(), "$GENE$");
	        		}
	        		else sentence = sentence.replace(mat2.group(), "$"+mat2.group()+"|GENE$");
	        	}
	        	
	        	String var = "";
	        	while(mat.find()) {//get variant
//	        		System.out.println(mat.group());
	        		if(!var.equals("")){
	        			var = var+","+mat.group();
	        		}
	        		else var = mat.group();
	        		sentence = sentence.replace(mat.group(), "$"+mat.group().replace(" ", "_")+"|MUT$");
	        	}
	        	
	        	if(!sentence.contains("MUT$")) {//get mutation
	        		while(mat1.find()) {
//		        		System.out.println(mat1.group());
		        		sentence = sentence.replace(mat1.group(), "$"+mat1.group().replace(" ", "_")+"|MUT$");
		        	}
	        	}
	        	
	        	while(mat3.find()) {//remove bracket
	        		if(!mat3.group().contains("$")&&!mat3.group().matches("NCCN|PARP|HBOC")) {
	        			sentence = sentence.replaceAll(mat3.group(), "");
	        		}
	        	}
	        	if(!mat3.find()) {
	        		sentence = sentence.replaceAll("\\(|\\)", "");
	        	}

	        	sentence = removeNumbers(sentence).replaceAll("\\s+", " ").trim();
	        	
	        	while(sentence.contains("|MUT$s")) {
	            	sentence = sentence.replace("|MUT$s", "|MUT$");
	        	}
	        	while(sentence.contains("|MUT$|MUT$")) {
	            	sentence = sentence.replace("|MUT$|MUT$", "|MUT$");
	        	}
	        	while(sentence.contains("|GENE$s")) {
	            	sentence = sentence.replace("|GENE$s", "|GENE$");
	        	}
	        	while(sentence.contains("|GENE$|GENE$")) {
	            	sentence = sentence.replace("|GENE$|GENE$", "|GENE$");
	        	}
	        	while(sentence.contains("$$")) {
	        		sentence = sentence.replace("$$", "$");
	        	}

	        	if(!mat3.find()) {
	        		sentence = sentence.replaceAll("\\(|\\)", "");
	        	}
	        	
	        	pat3 = Pattern.compile("\\b(?i)PARP(i)?\\b");
	        	mat3 = pat3.matcher(sentence);
	        	String mod = "";
	        	while(mat3.find()) {
	        		sentence = sentence.replace(mat3.group(), "$PARP$").replace("$$", "$");
	        		mod = "PARP";
	        	}
	        	pat3 = Pattern.compile("\\b(?i)germ(-|\\s)?line\\b");
	        	mat3 = pat3.matcher(sentence);
	        	while(mat3.find()) {
	        		sentence = sentence.replace(mat3.group(), "$germline$").replace("$$", "$");
	        		if(mod.equals("")) {
	        			mod = "germline";
	        		}
	        		else if(!mod.contains("germline")) {
	        			mod = mod+",germline";
	        		}
	        	}
	        	pat3 = Pattern.compile("\\b(?i)somatic\\b");
	        	mat3 = pat3.matcher(sentence);
	        	while(mat3.find()) {
	        		sentence = sentence.replace(mat3.group(), "$somatic$").replace("$$", "$");
	        		if(mod.equals("")) {
	        			mod = "somatic";
	        		}
	        		else if(!mod.contains("somatic")) {
	        			mod = mod+",somatic";
	        		}
	        	}
	        	
	        	String [] sentences = sentence.split("\\.");
	        	sentence = "";
	        	for(String s : sentences) {
	        		System.out.println(s);
	        		if(s.contains("$")) {
	        			if(sentence.equals("")) {
	        				sentence = s;
	        			}
	        			else sentence = sentence+"."+s;
		        	}
	        	}

	        	Matcher mat4 = pat4.matcher(sentence);
	        	Matcher mat5 = pat5.matcher(sentence);
	        	Matcher mat6 = pat6.matcher(sentence);
	        	Matcher mat7 = pat7.matcher(sentence);
	        	Matcher mat8 = pat8.matcher(sentence);
	        	Matcher mat9 = pat9.matcher(sentence);
	        	Matcher mat10 = pat10.matcher(sentence);
	        	System.out.println(sentence);
	        	
	        	String pos = "";
	        	while(mat9.find()){
	        		pos = "Y";
	        	}

	        	String output = ""; 
	        	if(sentence.toLowerCase().startsWith("it showed no ")||sentence.toLowerCase().startsWith("no ")){
	        		if(sentence.toLowerCase().contains("gene")||mat.find()){
	        			topic = "negative";
	    	        	pat = Pattern.compile("\\$([A-Za-z0-9\\.*>_]+)\\|MUT\\$");
	    	        	mat = pat.matcher(sentence);
	    	        	while(mat.find()) {
	    	        		sentence = sentence.replace(mat.group(), "MUT");
	    	        	}
	    	        	
	    	        	pat = Pattern.compile("\\$([A-Za-z0-9/]+)\\|GENE\\$");
	    	        	mat = pat.matcher(sentence);
	    	        	while(mat.find()) {
	    	        		sentence = sentence.replace(mat.group(), "BRCA");
	    	        	}
	    	        	sentence = sentence.replaceAll("[\\p{Punct}\\d]+", " ").replaceAll("\\s{2,}"," ");
	            		out.add(gene+"\t"+sentence+"\t"+topic+"\t"+mod+"\t"+var);
	            	}
	            	else {
	    	        	pat = Pattern.compile("\\$([A-Za-z0-9\\.*>_]+)\\|MUT\\$");
	    	        	mat = pat.matcher(sentence);
	    	        	while(mat.find()) {
	    	        		sentence = sentence.replace(mat.group(), "MUT");
	    	        	}
	    	        	
	    	        	pat = Pattern.compile("\\$([A-Za-z0-9/]+)\\|GENE\\$");
	    	        	mat = pat.matcher(sentence);
	    	        	while(mat.find()) {
	    	        		sentence = sentence.replace(mat.group(), "BRCA");
	    	        	}
	    	        	sentence = sentence.replaceAll("[\\p{Punct}\\d]+", " ").replaceAll("\\s{2,}"," ");
	            		out.add(gene+"\t"+sentence+"\t"+topic+"\t"+mod+"\t"+var);
	            	}
	            	continue;
        		}
	        	else if(sentence.contains("$"+gene+"|GENE$ mutant")){
	        		topic = "mutant";
		        	pat = Pattern.compile("\\$([A-Za-z0-9\\.*>_]+)\\|MUT\\$");
		        	mat = pat.matcher(sentence);
		        	while(mat.find()) {
		        		sentence = sentence.replace(mat.group(), "MUT");
		        	}
		        	
		        	pat = Pattern.compile("\\$([A-Za-z0-9/]+)\\|GENE\\$");
		        	mat = pat.matcher(sentence);
		        	while(mat.find()) {
		        		sentence = sentence.replace(mat.group(), "BRCA");
		        	}
		        	sentence = sentence.replaceAll("[\\p{Punct}\\d]+", " ").replaceAll("\\s{2,}"," ");
	        		out.add(gene+"\t"+sentence+"\t"+topic+"\t"+mod+"\t"+var);
	        		continue;
	        	}
	        	else if(sentence.toLowerCase().contains("insurance")||sentence.toLowerCase().contains("medicaid")||sentence.toLowerCase().contains("medicare")){
	        		topic = "insurance";
		        	pat = Pattern.compile("\\$([A-Za-z0-9\\.*>_]+)\\|MUT\\$");
		        	mat = pat.matcher(sentence);
		        	while(mat.find()) {
		        		sentence = sentence.replace(mat.group(), "MUT");
		        	}
		        	
		        	pat = Pattern.compile("\\$([A-Za-z0-9/]+)\\|GENE\\$");
		        	mat = pat.matcher(sentence);
		        	while(mat.find()) {
		        		sentence = sentence.replace(mat.group(), "BRCA");
		        	}
		        	sentence = sentence.replaceAll("[\\p{Punct}\\d]+", " ").replaceAll("\\s{2,}"," ");
	        		out.add(gene+"\t"+sentence+"\t"+topic+"\t"+mod+"\t"+var);
	        		continue;
	        	}
        		else if(mat4.find()){
        			topic = "order";
    	        	pat = Pattern.compile("\\$([A-Za-z0-9\\.*>_]+)\\|MUT\\$");
    	        	mat = pat.matcher(sentence);
    	        	while(mat.find()) {
    	        		sentence = sentence.replace(mat.group(), "MUT");
    	        	}
    	        	
    	        	pat = Pattern.compile("\\$([A-Za-z0-9/]+)\\|GENE\\$");
    	        	mat = pat.matcher(sentence);
    	        	while(mat.find()) {
    	        		sentence = sentence.replace(mat.group(), "BRCA");
    	        	}
    	        	sentence = sentence.replaceAll("[\\p{Punct}\\d]+", " ").replaceAll("\\s{2,}"," ");
	        		out.add(gene+"\t"+sentence+"\t"+topic+"\t"+mod+"\t"+var);
	        		continue;
	        	}
        		else if(((mat5.find()&&!sentence.contains("old"))||mat6.find())&&pos.equals("")){
        			topic = "discuss";
    	        	pat = Pattern.compile("\\$([A-Za-z0-9\\.*>_]+)\\|MUT\\$");
    	        	mat = pat.matcher(sentence);
    	        	while(mat.find()) {
    	        		sentence = sentence.replace(mat.group(), "MUT");
    	        	}
    	        	
    	        	pat = Pattern.compile("\\$([A-Za-z0-9/]+)\\|GENE\\$");
    	        	mat = pat.matcher(sentence);
    	        	while(mat.find()) {
    	        		sentence = sentence.replace(mat.group(), "BRCA");
    	        	}
    	        	sentence = sentence.replaceAll("[\\p{Punct}\\d]+", " ").replaceAll("\\s{2,}"," ");
	        		out.add(gene+"\t"+sentence+"\t"+topic+"\t"+mod+"\t"+var);
	        		continue;
        		}
        		else if((mat8.find()||(!sentence.toLowerCase().contains("cardio")&&mat7.find())&&pos.equals(""))){
	        		topic = "estimate";
		        	pat = Pattern.compile("\\$([A-Za-z0-9\\.*>_]+)\\|MUT\\$");
		        	mat = pat.matcher(sentence);
		        	while(mat.find()) {
		        		sentence = sentence.replace(mat.group(), "MUT");
		        	}
		        	
		        	pat = Pattern.compile("\\$([A-Za-z0-9/]+)\\|GENE\\$");
		        	mat = pat.matcher(sentence);
		        	while(mat.find()) {
		        		sentence = sentence.replace(mat.group(), "BRCA");
		        	}
		        	sentence = sentence.replaceAll("[\\p{Punct}\\d]+", " ").replaceAll("\\s{2,}"," ");
	        		out.add(gene+"\t"+sentence+"\t"+topic+"\t"+mod+"\t"+var);
	        		continue;
	        	}
        		else if(mat10.find()&&pos.equals("")){
	        		topic = "estimate+discuss";
		        	pat = Pattern.compile("\\$([A-Za-z0-9\\.*>_]+)\\|MUT\\$");
		        	mat = pat.matcher(sentence);
		        	while(mat.find()) {
		        		sentence = sentence.replace(mat.group(), "MUT");
		        	}
		        	
		        	pat = Pattern.compile("\\$([A-Za-z0-9/]+)\\|GENE\\$");
		        	mat = pat.matcher(sentence);
		        	while(mat.find()) {
		        		sentence = sentence.replace(mat.group(), "BRCA");
		        	}
		        	sentence = sentence.replaceAll("[\\p{Punct}\\d]+", " ").replaceAll("\\s{2,}"," ");
	        		out.add(gene+"\t"+sentence+"\t"+topic+"\t"+mod+"\t"+var);
	        		continue;
	        	}
	        	else {
	        		
   	        		String geneName = "";
		        	pat = Pattern.compile("\\$([A-Za-z0-9\\.*>_]+)\\|MUT\\$");
		        	mat = pat.matcher(sentence);

		        	pat1 = Pattern.compile("\\$([A-Za-z0-9/]+)\\|GENE\\$");
		        	mat1 = pat1.matcher(sentence);
		        	while(mat1.find()) {
		        		gene = mat1.group();
		        		geneName = mat1.group(1);
		        	}
		        	
		        	int min = 5;
	        		while(mat.find()){//gene - mutation relationship
	        			if(!getRelationship2(output, sentence, min, gene, mat.group()).equals("")) {
		        			output = getRelationship2(output, sentence, min, gene, mat.group());
		        			min = Integer.valueOf(output.split("\t")[0]);
	        			}
		        		while (countOccurences(sentence, gene)>1){
		        			sentence = sentence.replaceFirst(geneName, "").replaceAll("\\s+", " ");
		        			if(!getRelationship2(output, sentence, min, gene, mat.group()).equals("")&&Integer.valueOf(getRelationship2(output, sentence, min, gene, mat.group()).split("\t")[0])<min) {
			        			output = getRelationship2(output, sentence, min, gene, mat.group());
		        				min = Integer.valueOf(output.split("\t")[0]);
		        			}
		        		}
		        	}
	        	}
	        	
	        	String output1="";
	        	if(!output.equals("")) {
	        		System.out.println(output);
	        		int min = 5;
	        		mat9 = pat9.matcher(sentence);
		        	while(mat9.find()){
		        		if(!getRelationship3(output1, sentence, min, gene, output.split("\t")[2], mat9.group()).equals("")) {
		        			output1 = getRelationship3(output1, sentence, min, gene, output.split("\t")[2], mat9.group());
		        			min = Integer.valueOf(output1.split("\t")[0]);
		        		}
		        		while (countOccurences(sentence, mat9.group())>1){
		        			sentence = sentence.replaceFirst(mat9.group(), "").replaceAll("\\s+", " ");
		        			if(!getRelationship3(output1, sentence, min, gene, output.split("\t")[2], mat9.group()).equals("")&&Integer.valueOf(getRelationship3(output1, sentence, min, gene, output.split("\t")[2], mat9.group()).split("\t")[0])<min) {
			        			output1 = getRelationship3(output1, sentence, min, gene, output.split("\t")[2], mat9.group());
			        			min = Integer.valueOf(output1.split("\t")[0]);
			        		}
		        		}
		        	}
		        	
		        	if(!output1.equals("")) {
		        		
		        		System.out.println(output1);
		        		
		        		if(output1.split("\t")[1].toLowerCase().equals("vus")||output1.split("\t")[1].toLowerCase().startsWith("un")) {
		        			topic = "vus";
		        		}
		        		else if(output1.split("\t")[1].toLowerCase().equals("negative")||output1.split("\t")[1].toLowerCase().startsWith("no")||output1.split("\t")[1].toLowerCase().startsWith("not")) {
		        			topic = "negative";
		        		}
		        		else if(output1.split("\t")[1].toLowerCase().equals("would")||output1.split("\t")[1].toLowerCase().equals("if")||output1.split("\t")[1].toLowerCase().equals("might")){
		        			topic = "possible";
		        		}
		        		else if(output1.split("\t")[1].toLowerCase().equals("known")||output1.split("\t")[1].toLowerCase().equals("pathogenic")||output1.split("\t")[1].toLowerCase().equals("positive")){
		        			topic = "positive";
		        		}
		        		else {
		        			topic = "positive";
		        		}
	        		}
		        	else {
		        		topic = "positive";
		        	}
		        	System.out.println(topic);
	        	}
	        	else {
	        		
	        		int min = 5;
	        		mat9 = pat9.matcher(sentence);
		        	while(mat9.find()){
		        		System.out.println(mat9.group());
		        		if(!getRelationship2(output, sentence, min, gene, mat9.group()).equals("")) {
		        			output = getRelationship2(output, sentence, min, gene, mat9.group());
		        			min = Integer.valueOf(output.split("\t")[0]);
		    			}
		        		while (countOccurences(sentence, mat9.group())>1){
		        			sentence = sentence.replaceFirst(mat9.group(), "").replaceAll("\\s+", " ");
		        			if(!getRelationship2(output, sentence, min, gene, mat9.group()).equals("")) {
		            			output = getRelationship2(output, sentence, min, gene, mat9.group());
		            			min = Integer.valueOf(output.split("\t")[0]);
		        			}
		        		}
		        	}

		        	if(!output.equals("")) {
		        		
		        		System.out.println(output);
		        		
		        		if(output.split("\t")[2].toLowerCase().equals("positive")||output.split("\t")[2].toLowerCase().equals("pathogenic")) {
		        			topic = "positive";
		        		}
		        		else if(output.split("\t")[2].toLowerCase().equals("vus")) {
		        			topic = "vus";
		        		}
		        		else if(output.split("\t")[2].toLowerCase().equals("would")||output.split("\t")[2].toLowerCase().equals("if")||output.split("\t")[2].toLowerCase().equals("might")){
		        			topic = "possible";
		        		}
		        		else {
		        			topic = "negative";
		        		}
	        		}
		        	else {
		        		mat9 = pat9.matcher(sentence);
		        		System.out.println(countOccurences(sentence, "$GENE$"));
			        	while(mat9.find()) {
			        		if(countOccurences(sentence, "$GENE$")>=4&&(mat9.group().contains("no")||mat9.group().contains("negative"))) {
			        			topic = "negative";
			        		}
			        	}
		        	}
	        	}
	        	
//	        	gene = line.split("\t")[7];
	        	pat = Pattern.compile("\\$([A-Za-z0-9\\.*>_]+)\\|MUT\\$");
	        	mat = pat.matcher(sentence);
	        	while(mat.find()) {
	        		sentence = sentence.replace(mat.group(), "MUT");
	        	}
	        	
	        	pat = Pattern.compile("\\$([A-Za-z0-9/]+)\\|GENE\\$");
	        	mat = pat.matcher(sentence);
	        	while(mat.find()) {
	        		sentence = sentence.replace(mat.group(), "BRCA");
	        	}
	        	sentence = sentence.replaceAll("[\\p{Punct}\\d]+", " ").replaceAll("\\s{2,}"," ");
	        	out.add(gene.replace("$", "").replace("|GENE", "")+"\t"+sentence+"\t"+topic+"\t"+mod+"\t"+var);
//	        	System.out.println(gene+"\t"+sentence+"\t"+topic);
	        }
	    }
	}
    
	public static String removeStopWords(String input){
		
		ArrayList<String> wordsList = new ArrayList<String>();
	    input = input.trim().replaceAll("\\s+", " ");
        String[] words = input.split(" ");

        for (String word : words) {
            wordsList.add(word);
        }
        
        // remove stop words here from the temp list
        for (int i = 0; i < wordsList.size(); i++) {
            // get the item as string
            for (int j = 0; j < stopwords.length; j++) {
                if (wordsList.contains(stopwords[j])) {
                    wordsList.remove(stopwords[j]);//remove it
                }
            }
        }
        
        // output
        String output = "";
        if(wordsList.size()>0) {
    		output = wordsList.get(0);
            for (int i = 1; i < wordsList.size(); i++) {
            	output = output + " " + wordsList.get(i);
            }
        }
        return output;
	}
	
	public static String removeNumbers(String input){
		
		ArrayList<String> wordsList = new ArrayList<String>();
	    input = input.trim().replaceAll("\\s+", " ");
        String[] words = input.split(" ");

        for (String word : words) {
        	if(!word.matches("\\d+")&&!word.matches("[A-Za-z]")) {
        		wordsList.add(word);
        	}
        }
        
        // output
        String output = "";
        if(wordsList.size()>0) {
    		output = wordsList.get(0);
            for (int i = 1; i < wordsList.size(); i++) {
            	output = output + " " + wordsList.get(i);
            }
        }
        return output;
	}
	
	public static String getRelationship2(String output, String sentence, int min, String aa, String bb) 
	{
		int a = searchedWord(sentence, aa);
		int b = searchedWord(sentence, bb);
		if(Math.abs(a-b)<min&&b!=-1) {
			output = Math.abs(a-b)+"\t"+aa+"\t"+bb+"\t"+sentence;
			//min gene mut sentence
		}
		
		return output;
	}
	
	public static String getRelationship3(String output1, String sentence, int min, String aa, String bb, String cc) 
	{
		sentence = sentence.replace(cc, cc.replaceAll(" ", "_"));
		cc = cc.replaceAll(" ", "_");
		
		int a = searchedWord(sentence, aa);
		int b = searchedWord(sentence, bb);
		int c = searchedWord(sentence, cc);
		if((Math.abs(b-c)<min&&Math.abs(b-c)<=Math.abs(a-c))&&c!=-1) {
			min = b-c;
			output1 = Math.abs(b-c)+"\t"+cc+"\t"+aa+"\t"+bb;
		}
		else if((Math.abs(a-c)<min&&Math.abs(b-c)>Math.abs(a-c))&&c!=-1) {
			min = a-c;
			output1 = Math.abs(a-c)+"\t"+cc+"\t"+aa+"\t"+bb;
		}
		
		return output1;
	}
	
	public static int countOccurences(String str, String word) 
	{
	    // split the string by spaces in a
	    String a[] = str.split(" ");
	
	    // search for pattern in a
	    int count = 0;
	    for (int i = 0; i < a.length; i++) 
	    {
	    // if match found increase count
	    if (word.equals(a[i]))
	        count++;
	    }
	
	    return count;
	}
	
	public static int searchedWord(String sentence, String searchedWord) {
		
	    if (sentence == null || searchedWord == null) throw new IllegalArgumentException("May not be null");
	    String[] words = sentence.split(" ");
	    int searchedIndex = -1;
	    for (int i=0; i<words.length; i++) {
	        if (words[i].equals(searchedWord)) searchedIndex = i;
	        if (searchedIndex!=-1) break;
	    }

	    int output = searchedIndex;
//	    System.out.println(output);
	    return output;
	}
}
