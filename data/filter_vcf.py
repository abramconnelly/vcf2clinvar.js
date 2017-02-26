#!/usr/bin/python
#
# Create a new VCF from clinvar.vcf
# that only has the information fields
# we actually need.
#

with open("clinvar.vcf") as fp:
  for line in fp:
    if len(line)==0: continue
    if line[0]=='#':
      print line.strip()
      continue
    line = line.strip()
    fields = line.split("\t")
    chrom = fields[0]
    pos = fields[1]
    rsid = fields[2]
    ref = fields[3]
    alt = fields[4]
    x = fields[5]
    y = fields[6]
    info_fields = fields[7].split(";")

    m_field = {}

    #foi = { "CLNDBN":True, "CLNACC":True, "CAF":True }
    #afoi = [ "CLNDBN", "CLNACC", "CAF", ]

    foi = {'CAF':True,'CLNALLE':True, 'CLNDSDB':True, 'CLNDSDBID':True,
           'CLNACC':True, 'CLNDBN':True, 'CLNSIG':True,
           'CLNHGVS':True, 'CLNSRC':True, 'CLNSRCID':True }
    afoi = ['CAF','CLNALLE', 'CLNDSDB', 'CLNDSDBID', 'CLNACC', 'CLNDBN', 'CLNSIG', 'CLNHGVS', 'CLNSRC', 'CLNSRCID']

    for k in foi:
      m_field[k] = ""

    for ifield in info_fields:
      ifields = ifield.split("=")
      if ifields[0] in foi:
        m_field[ifields[0]] = ifields[1]

    info_a = []
    for a in afoi:
      if m_field[a] != "":
        info_a.append( a + "=" + m_field[a] )

    info_str = ";".join(info_a)

    #print "\t".join( [ chrom, pos, rsid, ref, alt, m_field[afoi[0]], m_field[afoi[1]], m_field[afoi[2]] ] )
    print "\t".join( [ chrom, pos, ".", ref, alt, x, y, info_str ] )

