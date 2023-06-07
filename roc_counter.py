# Script 4: ROC Curves ----------------------------------------------------------------

# Study Details -----------------------------------------------------------
# Title: "Assessing the performance of ten pathogenicity prediction algorithms on a dataset of missense brca2 and brca2 mutations."
# Author: Bertalan Takács
# Date: 16/05/2022
# Laboratory: HCEMM-BRC Mutagenesis and Carcinogenesis Research Group, Institute of Genetics, Biological Research Centre, Szeged H-6726, Hungary.

from sys import argv
import os
import numpy as np
from sklearn import metrics
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager

print(matplotlib.matplotlib_fname())

def collect_data_from_files(folder):
	file_list = os.listdir(folder)
	csvs = []
	predictors_results = {}
	for f in file_list:
		if f[-4:] == ".csv":
			csvs.append(f)
	for filename in csvs:
		with open(folder + filename,"r") as prediction:
			predictors_results[filename[:-10]] = {}
			predictors_results[filename[:-10]]["minimum"] = 0
			predictors_results[filename[:-10]]["maximum"] = 0
			predictors_results[filename[:-10]]["mutations"] = []
			for line in prediction:
				l = line.rstrip().split("\t")
				try:
					predictors_results[filename[:-10]]["mutations"].append([l[1],float(l[2])])
					if float(l[2]) < predictors_results[filename[:-10]]["minimum"]:
						predictors_results[filename[:-10]]["minimum"] = float(l[2])
					elif float(l[2]) > predictors_results[filename[:-10]]["maximum"]:
						predictors_results[filename[:-10]]["maximum"] = float(l[2])
				except:
					continue
	print(list(predictors_results))
	#Upper and lower bounds, what marks the pathogen?
	predictors_results["SIFT"]["limits"] = [0,1,"lower"]
	predictors_results["HumDiv"]["limits"] = [0,1,"higher"]
	predictors_results["HumVar"]["limits"] = [0,1,"higher"]
	predictors_results["SNPs_GO"]["limits"] = [0,1,"higher"]
	predictors_results["PROVEAN"]["limits"] = [-10.57,0,"higher"]
	predictors_results["PMut"]["limits"] = [0,1,"higher"]
	predictors_results["PhD-SNP"]["limits"] = [0,1,"higher"]
	predictors_results["META-SNP"]["limits"] = [0,1,"higher"]
	predictors_results["PredictSNP"]["limits"] = [0,1,"higher"]
	predictors_results["PANTHER"]["limits"] = [0,1,"higher"]
	return predictors_results
	
def calculate_scores(predictors,chunks):
	chunks = int(chunks)
	scores = {}
	for key in predictors:
		#TPR = TP/TP+FN
		#FPR = FP/FP+TN
		scores[key] = []
		try:
			chunk_sizes = abs(predictors[key]["limits"][0]-predictors[key]["limits"][1])/int(chunks)
		except:
			break
		if key not in ["SIFT","PROVEAN","SIFT-alternative"]:
			for i in range(1,chunks+1):
				tp = 0
				tn = 0
				fp = 0
				fn = 0
				threshold = i*chunk_sizes
				try:
					for l in predictors[key]["mutations"]:
						if l[1] > threshold:
							if l[0] == "pathogenic":
								tp += 1
							elif l[0] == "benign":
								fp += 1
						elif l[1] <= threshold:
							if l[0] == "pathogenic":
								fn += 1
							elif l[0] == "benign":
								tn += 1
					scores[key].append([1-(tn/(fp+tn)),(tp/(tp+fn))])
				except:
					continue
		elif key in ["SIFT","SIFT-alternative"]:
			chunk_sizes = abs(predictors[key]["maximum"]-predictors[key]["minimum"])/int(chunks)
			for i in range(-1,chunks+1):
				tp = 0
				tn = 0
				fp = 0
				fn = 0
				threshold = i*chunk_sizes
				for l in predictors[key]["mutations"]:
					if l[1] <= threshold:
						if l[0] == "pathogenic":
							tp += 1
						elif l[0] == "benign":
							fp += 1
					elif l[1] > threshold:
						if l[0] == "pathogenic":
							fn += 1
						elif l[0] == "benign":
							tn += 1
				scores[key].append([1-(tn/(fp+tn)),(tp/(tp+fn))])
		elif key == "PROVEAN":
			for i in range(1,chunks+1):
				tp = 0
				tn = 0
				fp = 0
				fn = 0
				threshold = i*chunk_sizes
				for l in predictors[key]["mutations"]:
					if l[1] <= -threshold:
						if l[0] == "pathogenic":
							tp += 1
						elif l[0] == "benign":
							fp += 1
					elif l[1] > -threshold:
						if l[0] == "pathogenic":
							fn += 1
						elif l[0] == "benign":
							tn += 1
				scores[key].append([1-(tn/(fp+tn)),(tp/(tp+fn))])
	return scores



"""
roc_auc = metrics.auc(fpr, tpr)

plt.figure()
plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.show()"""

def print_output(scores, folder):
	auc_scores = {}
	"""tpr_panther = []
	fpr_panther = []
	tpr_humdiv = []
	fpr_humdiv = []
	for eredm in scores["PANTHER"]:
			tpr_panther.append(eredm[1])
			fpr_panther.append(eredm[0])
	for eredm in scores["HumDiv"]:
			tpr_humdiv.append(eredm[1])
			fpr_humdiv.append(eredm[0])	
	tpr_panther.insert(0,1)
	fpr_panther.insert(0,1)
	roc_auc_panther = metrics.auc(np.array(tpr_panther),np.array(fpr_panther))
	tpr_humdiv.insert(0,1)
	fpr_humdiv.insert(0,1)
	roc_auc_humdiv = metrics.auc(np.array(tpr_humdiv),np.array(fpr_humdiv))
	plt.figure()
	plt.plot(fpr_panther, tpr_panther, label='Panther ROC curve (area = %0.2f)' % float(1-roc_auc_panther))
	plt.plot(fpr_humdiv, tpr_humdiv, label='HumDiv ROC curve (area = %0.2f)' % float(1-roc_auc_humdiv))
	plt.plot([0, 1], [0, 1], 'k--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic Panther and HumDiv')
	plt.legend(loc="lower right")
	plt.show()
	for key in scores:
		tpr = []
		fpr = []
		for eredm in scores[key]:
			tpr.append(eredm[1])
			fpr.append(eredm[0])
		tpr.sort(reverse = True)
		fpr.sort(reverse = True)
		tpr.insert(0,1)
		fpr.insert(0,1)
		tpr = np.array(tpr)
		fpr = np.array(fpr)
		roc_auc = metrics.auc(fpr, tpr)
		plt.figure()
		plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
		plt.plot([0, 1], [0, 1], 'k--')
		plt.xlim([0.0, 1.0])
		plt.ylim([0.0, 1.05])
		plt.xlabel('False Positive Rate')
		plt.ylabel('True Positive Rate')
		plt.title('Receiver operating characteristic ' + key)
		plt.legend(loc="lower right")
		plt.show()"""
	plt.rc('font', family='Helvetica')
	plt.figure()
	
	for key in scores:
		auc_scores[key] = {"tpr":[0,1],"fpr":[0,1], "auc" : 0}
		for eredm in scores[key]:
			auc_scores[key]["tpr"].append(eredm[1])
			auc_scores[key]["fpr"].append(eredm[0])
		auc_scores[key]["tpr"].sort(reverse = True)
		auc_scores[key]["fpr"].sort(reverse = True)
		auc_scores[key]["tpr"] = np.array(auc_scores[key]["tpr"])
		auc_scores[key]["fpr"] = np.array(auc_scores[key]["fpr"])
		auc_scores[key]["auc"] = metrics.auc(auc_scores[key]["tpr"],auc_scores[key]["fpr"])
		print(key + " AUC: " + str(1-auc_scores[key]["auc"]))
		try:
			plt.plot(auc_scores[key]["fpr"], auc_scores[key]["tpr"], label= key + ': %0.2f' % float(1-auc_scores[key]["auc"]))
		except:
			print(auc_scores[key])
			continue
		with open(folder + key + "_tpr_fpr.tsv","w") as out:
			for l in scores[key]:
				out.write("\t".join([str(x) for x in l]) + "\n")
	hfont = {'fontname':'Helvetica'}
	plt.plot([0, 1], [0, 1], 'k--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate', {'fontname':'Helvetica', "fontsize":11})
	plt.ylabel('True Positive Rate', {'fontname':'Helvetica', "fontsize":11})
	plt.title('ROC - BRCA', {'fontname':'Helvetica', "fontsize":13})
	plt.legend(loc="lower right")
	#plt.show()
	
	plt.rcParams['font.size'] = 12
	plt.savefig(folder + folder[0:-1] + ".png", dpi=500)

print_output(calculate_scores(collect_data_from_files(argv[1]),argv[2]),argv[1])
