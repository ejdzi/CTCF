def clusters(): #tworzenie pliku z wybranymi Raw PETsami, ktore lapia sie w obszary kotwic 
	f = open("reprezentacja_anchors.bed","r+") #"GM12878_anchors.txt","r+"
	g = open("raw_PETs_anchors.txt","r+") #"GSM1872886_GM12878_CTCF_PET_clusters.txt","r+"
	slownik={}
	h = open("raw_PETs_reprezentacja.txt","w+") #ten plik to Raw PETs "raw_PETs_anchors.txt","w+"
	for line in f:
		l = line.split()
		if l[0] not in slownik.keys():
			slownik[l[0]] = [(l[1],l[2])]
		elif l[0] in slownik.keys():
			slownik[l[0]].append((l[1],l[2]))
	for line in g:
		zmienna_p = False
		zmienna_k = False
		kotwica_p = None
		kotwica_k = None
		l = line.split()
		if l[0] in slownik.keys():
			lista = slownik[l[0]]
			for i in lista:
				if int(l[1]) >= int(i[0]) and int(l[2]) <=int(i[1]):
					kotwica_p = i
					zmienna_p = True
				if int(l[4]) >= int(i[0]) and int(l[5]) <= int(i[1]):
					kotwica_k = i
					zmienna_k = True
			if zmienna_p and zmienna_k:
				#h.write(l[0] + "\t" + l[1] + "\t" + l[2] + "\t" + l[4] + "\t" + l[5] +"\t" + l[6] +"\t" + str(kotwica_p) + "\t" +  str(kotwica_k) + "\n")
				h.write(l[0] + "\t" + l[1] + "\t" + l[2] + "\t" + l[3] + "\t" + l[4] +"\t" + l[5] + "\n") #"\t" + str(kotwica_p) + "\t" +  str(kotwica_k) + "\n")
	return slownik


def ctcf(): #tworzenie slownika wszystkich znalezionych CTCFow w obszarach kotwic --> dla reprezentacji sekwencji kotwic!
	f = open("reprezentacja_ctcf.bed","r+") # to jest dla reprezentacji
	slownik = {}
	for line in f:
		l = line.split()
		ide = l[0] + "\t" + l[1] + "\t" + l[2] 
		if ide not in slownik.keys():
			slownik[ide] = [(l[3],l[4],l[5])]
		elif ide in slownik.keys():
			slownik[ide].append((l[3],l[4],l[5]))
	return slownik

#trzeba teraz sprawdzac dla kazdej linii, czy w zakresie ich kotwic obu naraz znajduja sie jakis CTCF, jesli nie, to wtedy problem, niech sie takie wypsiza. 
# a jesli tak, no to ekstra --. to jest do zrobienia





def dynamic(l=[(0,5,15),(4,9,18),(8,21,19),(10,20,12),(25,30,25)]): #posortowane po poczatkach!
	lista = l[:] 
	scores = []
	for i in lista:	
		scores.append(i[2])
	s = [0] * len(scores)
	s[len(scores)-1] = scores[-1]
	indeksy=[len(scores)-1]
	wyniki = [0] * (len(lista))
	wyniki[len(lista)-1]=len(lista)-1
	for i in range(len(scores)-2,-1,-1):
		if lista[i][1] < lista[i+1][0]:
			s[i] = max(scores[i]+s[i+1],s[i+1])
			indeksy.append(i)
			wyniki[i]=(i,i+1)
		else:	
			zmienna = False
			for k in range(len(indeksy)-1,-1,-1):
				if lista[indeksy[k]][0] > lista[i][1]:
					j = indeksy[k]
					s[i] = max(scores[i] + s[j], s[j])
					zmienna = True
					wyniki[i] = (i,j)
					break
				if not zmienna:
					s[i] = scores[i]
					wyniki[i]=i
			indeksy.append(i)
	maks = 0
	i = -1
	for i in range(len(s)):
		if s[i] > maks:
			maks = s[i]
			m = i
	interwaly_wybrane = lokal(m,wyniki,lista=[])
	interwaly=[]
	for i in interwaly_wybrane:
		interwaly.append(lista[i])
	return interwaly #, s
		


def lokal(m,w=[], lista=[]): #daje indeksy poprawnych interwalow
	if type(w[m]) != tuple:
		lista.append(w[m])
		return lista
	else:
		lista.append(w[m][0])
		return lokal(w[m][1],w,lista) 

	
def niezachodzace(slownik):
	s = {}
	w = open("niezachodzace_ctcf_reprezentacja.bed","w+")
	for k in slownik.keys():
		lista = slownik[k]
		l = sorted(lista,key=lambda x: x[0])
		ctcfy= dynamic(l)
		kk = k.split("\t")[0]
		if kk not in s.keys():
			s[kk] = ctcfy
		elif kk in s.keys():
			s[kk].append(ctcfy)
	return s

'''teraz jak mam slownik z niezachodzacymi ctcfami dla kazdej kotwicy, musze przejsc po pliku z petsami od reprezentacji i sprawdzic, czy dla kazdej linii 
koniec lub poczatek petsa lapie sie do ktoregos zakresu ctcfa'''



def pets_vs_ctcf():
	slownik = niezachodzace(ctcf())
	f = open("raw_PETs_reprezentacja.txt","r+")
	for line in f:
		l = line.split()
		if l[













