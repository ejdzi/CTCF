import Bio
from Bio import motifs
from Bio import Seq

def gdzie_ctcf():
	f = open("GM12878reprezentacja_bed.faa","r+")
	print (f)
	g = open("reprezentacja_ctcf.bed","w+")
	motyw = motifs.read(open("ctcf_jaspar.txt"),"jaspar")
	for line in f:
		if line[0]==">":
			idd = line[1:-1]
			lista_idd=idd.split()
			print lista_idd
		l = f.next()
		sekwencja = l[:-1]
		s = Seq.Seq(sekwencja)
		s.alphabet = motyw.alphabet
		#s += s.reverse_complement()
		for x in motyw.pssm.search(s,threshold = 0):
			#print x
			if int(x[0]) < 0:
				start = int(lista_idd[2]) + int(x[0]) -1
				stop = start + 19 -1
				strand = "R"
			if int(x[0]) > 0:
				start = int(lista_idd[1]) + int(x[0]) -1
				stop = start + 19 - 1
				strand = "L"
			g.write(lista_idd[0] + "\t" + lista_idd[1] + "\t" + lista_idd[2] + "\t" + str(start) + "\t" + str(stop) +"\t" +  str(x[1]) + "\t" + strand + "\n")
	return


def dlugosc(plik = 'GM12878sekwencje_bed.faa'):
	lista = []
	lista_sekwencji=[]
	lista_10000_20000 = []
	lista_1000_2500 = []
	lista_2500_5000 = []
	lista_5000_7500 = []
	f = open(plik,"r")
	mini = ('id',100000000000000000)
	maxi = ('id',0)
	int_150_1000 = 0
	int_1000_2000 = 0
	int_2000_3000 = 0
	int_3000_4000 = 0
	int_1000_2500 = 0
	int_2500_5000 = 0
	int_5000_7500 = 0
	int_7500_10000 = 0
	average = 0
	for line in f:
		if line[0] == ">":
			idd = line[1:-1]
		l=f.next()
		s = l[:-1]
		lista.append((idd,len(s)))
		if len(s) < mini[1]:
			mini = (idd,len(s))
		if len(s) > maxi[1]:
			maxi = (idd,len(s))
		if 150 < len(s) < 10000:
			int_150_1000 += 1
		if 10001 < len(s) < 20000:
			int_1000_2000 +=1
			lista_10000_20000.append(idd)
		if 20001 < len(s) < 30000:
			int_2000_3000 +=1
		if 30001 < len(s) < 40000:
			int_3000_4000 += 1 
		if 1000 < len(s) < 2500:
			int_1000_2500 +=1
			lista_1000_2500.append(idd)
		if 2501 < len(s) < 5000:
			int_2500_5000 +=1
			lista_2500_5000.append(idd)
		if 5001 < len(s) < 7500:
			int_5000_7500 += 1
			lista_5000_7500.append(idd)
		if 7501 < len(s) < 1000: 
			int_7500_10000 += 1
		average += len(s)
	lista_sekwencji.append(mini[0])
	lista_sekwencji.append(maxi[0])
	#print lista_sekwencji
	l1=(lista_10000_20000[:25])
	l2=(lista_1000_2500[:25])
	l3=(lista_2500_5000[:25])
	l4=(lista_5000_7500[:25])
	lista_sekwencji += l1+l2+l3+l4
	average = float(average / len(lista))
	#return mini, maxi, average, int_1000_2500, int_2500_5000, int_5000_7500, int_7500_10000
	#return len(lista_10000_20000), len(lista_1000_2500), len(lista_2500_5000), len(lista_5000_7500)
	return len(lista_sekwencji), lista_sekwencji # to mi daje reprezentantow do szukania ctcfow w nich


def reprezentacja_sekwencji():
	lista = dlugosc()[1]
	f = open("GM12878sekwencje_bed.faa","r+")
	g = open("GM12878reprezentacja_bed.faa","w+")
	for line in f:
		if line[0] == ">":
			idd = line[1:-1]
		if idd in lista:
			l=f.next()
			s = l[:-1]
			g.write(">"+idd + "\n")
			g.write(s + "\n")
	f.close()
	g.close()

#reprezentacja_ctcf_t0.01.bed
def zlicz():
	f = open("reprezentacja_ctcf_t0.01.bed","r+")
	lista = dlugosc()[1]
	lista_reprezentacja=[]
	slownik={}
	for l in lista:
		slownik[l]=0
	count = 0
	print len(slownik)
	for line in f:
		l = line.split()
		idd = l[0] + " " + l[1] + " " + l[2]
		if idd in slownik.keys():
			slownik[idd]+=1
	return slownik, len(slownik) #stad slownik gdzie idzie
	


def plik_dlugosci_ilosci():
	f = open("GM12878reprezentacja_bed.faa","r+")
	g = open("GM12878reprezentacja_zliczenia_t0.01.csv","w+")
	slownik = zlicz()[0]
	g.write("Chromosom i zakres" + "\t" + "Dlugosc sekwencji" + "\t" + "Ilosc CTCF" + "\n")
	for line in f:
		if line[0] == ">":
			idd = line[1:-1]
		l=f.next()
		s = l[:-1]
		g.write(idd + "\t" + str(len(s)) + "\t" + str(slownik[idd]) + "\n")
	f.close()
	g.close()
		

def slownik_lokalizacji():
	f = open("reprezentacja_ctcf_t0.01.bed","r+")
	count = 0 
	slownik = {}
	lista = dlugosc()[1]
	for l in lista:
		slownik[l] = []
	for line in f:
		l = line.split()
		idd = l[0] + " " + l[1] + " " + l[2]
		if idd not in slownik.keys():
			print idd
		else:
			slownik[idd].append((l[3],l[4]))
	return slownik,len(slownik)



def zachodzi_napraw(slownik):
	zmienna = True
	for s in slownik.keys():
		lista = slownik[s]
		i = 0
		while i < (len(lista)-1):
			#print lista[i][0]
			if int(lista[i][1]) < int(lista[i+1][0]):
				i += 1
			else:
				zmienna = False
				start = lista[i][0]
				stop = lista[i+1][1]
				#print i, lista[i], lista[i+1]
				lista[i] = (start, stop)
				#print lista[i]
				lista=lista[:i+1]+lista[i+2:]
				#print lista
		slownik[s] = lista
				#print i, lista[:i+1], lista[i+2:]
				
	return slownik, zmienna, len(slownik)

def slownik_2_bed(slownik):
	w = open("reprezentacja_CTCF_niezachodzace_t0.01.bed","w+")
	for s in slownik.keys():
		lista = slownik[s]
		#print lista
		for i in lista:
			w.write(s + " " + i[0] + " " + i[1] + "\n")
			print i
	w.close()


def cluster():
	lista=[]
	f = open("reprezentacja_CTCF_niezachodzace_sorted.bed","r+")
	g = open("GM12878CTCF_clusters.txt","r+")
	h = open("reprezentacja_CTCF_clusters.bed","w+")
	for line in f:
		l = line.split()
		idd = l[0] + " " + l[1] + " " + l[2]
		if idd not in lista:
			lista.append(idd)
		#else:
		#	print len(lista)
	for line in g:
		l = line.split()
		idd_1 = l[0] + " " + l[1] + " " + l[2] 
		idd_2 = l[3] + " " + l[4] + " " + l[5]
		#print idd_1, idd_2
		if (idd_1 in lista):
			if idd_2 in lista:
				h.write(idd_1 + " " + idd_2 + " " + l[6] + "\n")
	return len(lista)
	
	

def dziury(slownik):
	w = open("CTCF_dziury_t0.01.bed","w+")
	for s in slownik:
		ss=s.split()
		dziury=[]
		lista=slownik[s]
		i = 0
		while i < (len(lista)-1):
			start = int(lista[i][1]) + 1
			stop = int(lista[i+1][0]) - 1
			dziury.append((str(start), str(stop)))
			w.write(ss[0] + " " + str(start) + " " + str(stop) + "\n")
			i+=1
	return



#ew zmienic okienko przeszukiwania
#moze lepiej znajdowac dziury miedzy CTCFami w kotwicach?
#mozna tez sprawdzic czy zakresy petclusters wypadaja w dziurach / ctcfach 

