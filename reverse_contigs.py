# coding: utf-8
contigs = open("contigs.fa")
contigs_rev = open("contigs_rev.fa", "w")
slownik = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
for linia in contigs:
    if linia.startswith(">"):
        contigs_rev.write(linia)
    else:
        linia = list(linia.strip())
        linia.reverse()
        for znak in linia:
            nowy_znak = slownik[znak]
            contigs_rev.write(nowy_znak)
        contigs_rev.write('\n')
