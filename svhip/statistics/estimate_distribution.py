
import testset_creation
import matplotlib.pyplot as plt
import sys
import os
import numpy
from scipy.stats import norm
from scipy.optimize import curve_fit
from astropy import modeling

dist_tuple = tuple()
argx = sys.argv

#Parameters:######
simulate =  False
samples = [80, 120, 160]
#smpl = 40
simsmples = [1]
mute = False
##################
'''
Creates comparision between alignment block distance distribution for raw data
and Rfam alignments shuffled by SISSIz before block creation. Arg1 and Arg2
therefore have to correspond.
Positive FIRST, negative SECOND. Name THIRD!
'''
def main():
	align_dict1 = {}
	align_dict2 = {}
	
	if len(argx) < 3:
		vector_dist, vector_n = read_points(str(argx[1]))
		mean_distances1, mean_distances2 = vector_dist[0], vector_dist[1]
		n_sequences1, n_sequences2 = vector_n[0], vector_n[1]
	else:
		setname = str(argx[3])
	
		for entry in os.listdir(argx[1]):
			filepath = os.path.join(argx[1], entry)
			align_dict1[str(entry)] = filepath
		for entry in os.listdir(argx[2]):
			filepath = os.path.join(argx[2], entry)
			align_dict2[str(entry)] = filepath
		
	for simsmpl in simsmples:
		for smpl in samples:
			#Determine and arrange values:	
			
			if not len(argx) < 3:
				mean_distances1, n_sequences1 = testset_creation.empirical_p_value(align_dict1, simulate, smpl, simsmpl, mute, rvalue =True)
				mean_distances2, n_sequences2 = testset_creation.empirical_p_value(align_dict2, simulate, smpl, simsmpl, mute, rvalue =True)
				write_data(mean_distances1, mean_distances2, n_sequences1, n_sequences2, str(setname) + '_' + str(smpl))
			
			print(n_sequences1)
			print(n_sequences2)
			
			distances_arranged1 = rearrange_gaussian(sorted(mean_distances1))
			distances_arranged2 = rearrange_gaussian(sorted(mean_distances2))
			#x_values = numpy.arange(len(mean_distances))

			#Plot measured values:
			fig, ax = plt.subplots()
			plt.title('Observed set: ' + str(setname) + '\n Distribution of averaged structural distances based on \n ' + str(smpl) +' alignment windows drawn from a seed alignment of 87 sequences. \n ')
			bins = numpy.arange(0, 100, 4)
			plt.xticks(numpy.arange(0, 100, 4))
			bin_heights1, bin_borders1, _ = plt.hist(mean_distances1, bins, color = 'g', label = 'positive set', alpha = 0.5)
			bin_heights2, bin_borders2, _ = plt.hist(mean_distances2, bins, color = 'r', label = 'randomized set', alpha = 0.5)

			'''
			#Try to fit curve for positive set
			mu1, std1 = norm.fit(distances_arranged1)
			mu2, std2 = norm.fit(distances_arranged2)
			
			try:
				bin_centers1 = bin_borders1[:-1] + numpy.diff(bin_borders1) / 2
				popt, _ = curve_fit(gauss, bin_centers1, bin_heights1, p0=[1., mu1, std1, 0.0])

				x_interval_for_fit = numpy.linspace(bin_borders1[0], bin_borders1[-1], 10000)
				plt.plot(x_interval_for_fit, gauss(x_interval_for_fit, *popt), label='positive fit', color = 'b')
				
			except Exception as e:
				print("Caught it. Nothing to see here.")
			
			#Try to fit curve for random set
			try:
				bin_centers2 = bin_borders2[:-1] + numpy.diff(bin_borders2) / 2
				popt, _ = curve_fit(gauss, bin_centers2, bin_heights2, p0=[1., mu2, std2, 0.0])

				x_interval_for_fit = numpy.linspace(bin_borders2[0], bin_borders2[-1], 10000)
				plt.plot(x_interval_for_fit, gauss(x_interval_for_fit, *popt), label='negative fit', color = 'b')
				
			except Exception as e:
				print("Caught it. Nothing to see here.")
			'''
			plt.legend()
			plt.grid()
			'''
			#Estimate a gaussian curve based on given data:
			#Currently 
			try:
				popt, pcov = curve_fit(gauss,x_values,ydata=distances_arranged, p0=[1,mu,std, 0.0])
				plt.plot(x_values,gauss(x_values,*popt),color ='r',label='fit')
			except Exception as e:
				print("There was a problem with curve fitting subroutine.")
			'''
			#Finish:
			plt.xlabel('avg. tree edit distance [single edit operations]')
			plt.ylabel('number of alignments in range')
			for tick in ax.get_xticklabels():
				tick.set_rotation(45)

			#Save the figure
			plt.savefig(str(setname) + '_smpl_'+ str(smpl)+'.png', format = 'png',  dpi = 1200)
			
			#Plot distances with respect to sequence number:
			plt.close()
			plt.title('Observed set: ' + str(setname) + '\n Average tree edit distance grouped by \n sequence count in alignment block. \n Sample size: ' + str(smpl))
			plt.plot(n_sequences1, mean_distances1, 'go', label = 'positive set', alpha = 0.5)
			plt.plot(n_sequences2, mean_distances2, 'ro', label = 'negative set', alpha = 0.5)
			plt.xlabel('number of sequences in alignment block')
			plt.ylabel('avg. tree edit distance [single edit operations]')
			plt.xticks(numpy.arange(0, 16, 1))
			for tick in ax.get_xticklabels():
				tick.set_rotation(45)
			plt.savefig(str(setname) + '_smpl_'+ str(smpl)+ '_grouped_by_seqs.png', format = 'png', dpi = 1200)

def rearrange_gaussian(data):
	data.reverse()
	sorted_data = []
	a = 0
	
	for point in data:
		a +=1 
		if a % 2 == 0:
			sorted_data.insert(0, point)
		else:
			sorted_data.append(point)

	return sorted_data

# Create a function which returns a Gaussian (normal) distribution.
def gauss(x, a, mu, sig, offset):
	y = a*numpy.exp(-numpy.power(x - mu, 2.) / (2 * numpy.power(sig, 2.))) + offset
	return y
	
def write_data(md1, md2, n1, n2, name):
	with open(str(name) + '.points', 'w+') as outf:
		for p in range(0, len(md1)):
			outf.write(str(md1[p]) + ':' +str(n1[p]) + ';')
		outf.write('\n')
		for p in range(0, len(md2)):
			outf.write(str(md2[p]) + ':' +str(n2[p]) + ';')
	print("Data points written to: " + str(name) + '.points') 
	return None
	
def read_points(filename):
	vector_dist = []
	vector_n = []
	with open(filename, 'r') as pf:
		for line in pf.readlines():
			ret_dist = []
			ret_n = []
			points = line.split(';')
			for e in points:
				e_split = e.split(':')
				ret_dist.append(e_split[0])
				ret_n.append(e_split[1])
			vector_dist.append(ret_dist)
			vector_n.append(ret_n)
	return vector_dist, vector_n

if __name__ == "__main__":
    main()
