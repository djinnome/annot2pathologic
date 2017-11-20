#!/usr/bin/env python
import argparse
from Bio import SeqIO
import os
import pandas as pd
import re
import operator
import numpy as np
import gffutils

class Annot2Pathologic:
    def __init__( self, annot):
        self.annot = annot
    def getProductType( self, geneId ):
        return 'P'
    def getName( self, geneId ):
        return self.annot.loc[geneId, 'Name']
    def getFunctions( self, geneId ):
        pass
    def getDBlinks( self, geneId ):
        pass
    
class PathologicEntry:
    attributes = ['ID','NAME','STARTBASE','ENDBASE',
                  'PRODUCT-TYPE','SYNONYM',
                  'DBLINK','GENE-COMMENT','GO',
                  'PRODUCT-ID','INTRON','functions', 'genetic_element']
    function_attributes = ['FUNCTION','EC','FUNCTION-COMMENT','FUNCTION-SYNONYM']
    repeat_attributes = ['SYNONYM', 'DBLINK','INTRON', 'GO']
    required_attributes = ['ID','NAME','STARTBASE','ENDBASE','functions','PRODUCT-TYPE']
    def __init__(self, **entry ):
        for att in entry:
            if att not in self.attributes:
                raise AttributeError('{} not a valid attribute'.format(att))
        for req in self.required_attributes:
            if req not in entry:
                raise AttributeError('Required attribute {} is missing from entry.'.format(req) )
        for funatt in entry['functions']:
            if 'FUNCTION' not in funatt:
                raise AttributeError('Required function attribute {} is missing from entry.'.format(funatt))
        self.entry = entry
    def __str__( self ):
        out = []
        for att in self.attributes:
            if att in self.entry:
                if att == 'functions':
                    for function_entry in  self.entry[att]:
                        for funatt in self.function_attributes:
                            if funatt in function_entry:
                                out.append('{}\t{}'.format(funatt, function_entry[funatt]))
                elif att in self.repeat_attributes and type(self.entry[att]) is list:
                    for repeated_att in self.entry[att]:
                        out.append('{}\t{}'.format(att,repeated_att))
                elif att == 'genetic_element':
                    pass
                else:
                    out.append('{}\t{}'.format(att, self.entry[att]))
        out.append('//\n')
        return '\n'.join(out)

class GFF2Pathologic:
    def __init__(self, gtf_file, dbfn='genome.db', force=True, merge_strategy='create_unique', **kwargs ):
        self.db = gffutils.create_db(gtf_file, dbfn=dbfn, force=force, merge_strategy=merge_strategy, **kwargs)
    def getGenes( self ):
        return self.db.features_of_type( 'gene', order_by='start')
    def getStartBase( self, gene ):
        if gene.strand == '+':
            return gene.start
        else:
            return gene.stop
    def get_non_cds( self, gene ):
        """This is actually only getting the non-cds part of the gene, not the introns"""
        introns = []
        cdss = []
        for cds in self.db.children( gene, featuretype='cds', order_by='start'):
            cdss.append( cds )
        
        if len(cdss) > 0 and gene.start < cdss[0].start:
            introns.append((gene.start, cdss[0].start -1))
        if len(cdss) > 1:
            for i in range(len(cdss) -1):
                introns.append((cdss[i].stop + 1, cdss[i+1].start - 1))
        if len(cdss) > 0 and gene.stop > cdss[-1].stop:
            introns.append((cdss[-1].stop + 1, gene.stop))
        return ['{:d}-{:d}'.format(startbase, endbase) for startbase, endbase in introns]
        
    def getEndBase( self, gene ):
        if gene.strand == '+':
            return gene.stop
        else:
            return gene.start
    def getId( self, gene ):
        return gene['portal_id'][0] + '_' + gene['proteinId'][0]

    def getGeneticElement( self, gene ):
        return gene.seqid

class JGIAnnot(Annot2Pathologic):
    def __init__( self, annot , sep='|'):
        self.annot = annot
        self.sep= sep
    def notEmpty( self, geneId, col ):
        return type(self.annot.loc[geneId, col]) is str
    def getEntry( self, geneId, col):
        return sorted( set([c for c in self.annot.loc[geneId, col].split(self.sep)]))
    def getName( self, geneId, names=['ECdef','kogdef', 'GOdef']):
        if geneId not in self.annot.index:
            return geneId
        else:
            for name in names:
                if self.notEmpty( geneId, name):
                    return self.getEntry(geneId, name)[0]
        return geneId
    def getFunctions( self, geneId ):
        functions = []
        if geneId not in self.annot.index:
            functions.append( dict(FUNCTION='ORF') )
        elif self.notEmpty(geneId, 'EC') and self.notEmpty(geneId, 'ECdef'):
            for ec, ecdef in zip(self.getEntry(geneId, 'EC'), self.getEntry(geneId, 'ECdef')):
                functions.append( dict(FUNCTION=ecdef, EC=ec))
        elif self.notEmpty( geneId, 'kogdef'):
            for kogdef in self.getEntry(geneId, 'kogdef'):
                functions.append( dict(FUNCTION=kogdef) )
        elif self.notEmpty( geneId, 'GOdef'):
            for godef in self.getEntry( geneId, 'GOdef'):
                functions.append( dict(FUNCTION=godef))
        return functions
    def getDBlinks( self, geneId ):
        dblinks = []
        if geneId not in self.annot.index:
            return dblinks
        if self.notEmpty( geneId, 'kog'):
            for kog in self.getEntry( geneId, 'kog'):
                dblinks.append('KOG:{}'.format(kog))
        if self.notEmpty( geneId, 'GO'):
            for go in self.getEntry( geneId, 'GO'):
                dblinks.append( go )
        return dblinks


class GTF2Pathologic:
    def __init__(self, gtf_file, dbfn='genome.db', force=True, merge_strategy='create_unique', **kwargs ):
        self.db = gffutils.create_db(gtf_file, dbfn=dbfn, force=force, merge_strategy=merge_strategy, **kwargs)
    def gene2entry( self, gene, geneKey  ):
        pe = {}
        functions=[]
        pe['ID'] = self.getId( gene )
        pe['genetic_element'] = self.getGeneticElement( gene )
        pe['STARTBASE'] = self.getStartBase( gene )
        pe['ENDBASE'] = self.getEndBase( gene )
        pe['INTRON'] = self.get_non_cds( gene )                
        pe['PRODUCT-TYPE'] = self.getProductType( gene )
        pe['NAME'] = pe['ID']
        pe['functions'] = self.getFunctions( gene )
        return pe
    def getProductType( self, gene ):
        return 'P'

    def get_entries( self, geneKey ):
        pathologic_entries = {}
        for gene in self.getGenes():
            pe = self.gene2entry( gene, geneKey )
            if pe['genetic_element'] in pathologic_entries:
                pathologic_entries[pe['genetic_element']].append( pe )
            else:
                pathologic_entries[pe['genetic_element']] = [pe]
        return pathologic_entries
    def generate_pathologic_files( self, pathologic_entries, pathologic_files={},template_file='genetic_element_{}.pf' , outputdir='.'):
        for genetic_element in pathologic_entries:
            if genetic_element not in pathologic_files:
                pathologic_files[genetic_element] = template_file.format(genetic_element )
                with open(os.path.join(outputdir,pathologic_files[genetic_element]), 'w') as pf:
                    for pe in sorted(pathologic_entries[genetic_element], key=operator.itemgetter('ID')):
                        pf.write(str(PathologicEntry(**pe)))
            else:
                with open(os.path.join(outputdir,pathologic_files[pe['genetic_element']]),'a') as pf:
                    for pe in sorted(pathologic_entries[genetic_element], key=operator.itemgetter('ID')):
                        pf.write(str(PathologicEntry(**pe)))
        return pathologic_files

    def generate_genetic_elements_file( self, pathologic_files,seq_files={}, circular={},outputdir='.'):
        with open(os.path.join(outputdir,'genetic-elements.dat'),'w') as out:
            for genetic_element in sorted(pathologic_files):
                out.write('ID\t{}\n'.format(genetic_element.replace('_','')))
                out.write('NAME\t{}\n'.format( genetic_element))
                out.write('TYPE\t{}\n'.format( ':CHRSM'))
                if genetic_element in circular:
                    out.write('CIRCULAR?\t{}\n'.format(circular[genetic_element]))
                else:
                    out.write('CIRCULAR\tN\n')
                out.write('ANNOT-FILE\t{}\n'.format(pathologic_files[genetic_element]))
                if genetic_element in seq_files:
                    out.write('SEQ-FILE\t{}\n'.format(seq_files[genetic_element]))
                out.write('//\n')

    def getGenes( self ):
        return self.db.features_of_type( 'transcript', order_by='start')
    def get_non_cds( self, gene ):
        """This is actually only getting the non-cds part of the gene, not the introns"""
        introns = []
        cdss = []
        for cds in self.db.children( gene, featuretype='cds', order_by='start'):
            cdss.append( cds )
        
        if len(cdss) > 0 and gene.start < cdss[0].start:
            introns.append((gene.start, cdss[0].start -1))
        if len(cdss) > 1:
            for i in range(len(cdss) -1):
                introns.append((cdss[i].stop + 1, cdss[i+1].start - 1))
        if len(cdss) > 0 and gene.stop > cdss[-1].stop:
            introns.append((cdss[-1].stop + 1, gene.stop))
        return ['{:d}-{:d}'.format(startbase, endbase) for startbase, endbase in introns]
    
    def get_introns_only_not_UTR( self, gene ):
        introns = []
        for intron in self.db.children( gene, featuretype='intron', order_by='start'):
            introns.append('{:d}-{:d}'.format(intron.start, intron.stop))
        return introns
        #exons = []
        #for exon in self.db.children( gene, featuretype='CDS', order_by='start'):
        #    exons.append((exon.start, exon.stop))
        #if len(exons) > 1:
        #    for i in range(len(exons) -1):
        #        introns.append((exons[i][1] + 1, exons[i+1][0] - 1))
        #return ['{:d}-{:d}'.format(startbase, endbase) for startbase, endbase in introns]
    def getStartBase( self, gene ):
        if gene.strand == '+':
            return gene.start
        else:
            return gene.stop
    def getEndBase( self, gene ):
        if gene.strand == '+':
            return gene.stop
        else:
            return gene.start
    def getId( self, gene ):
        return gene['transcript_id'][0]

    def getGeneticElement( self, gene ):
        return gene.seqid
    def getFunctions( self, geneId ):
        return [dict(FUNCTION='ORF')]
    
    
class GFFandAnnot2Pathologic:
    def __init__(self, gff, annot ):
        self.annot = annot
        self.gff = gff
    
    def gene2entry( self, gene, geneKey  ):
        pe = {}
        geneId = int(gene[geneKey][0])
        functions=[]
        pe['ID'] = self.gff.getId( gene )
        pe['genetic_element'] = self.gff.getGeneticElement( gene )
        pe['STARTBASE'] = self.gff.getStartBase( gene )
        pe['ENDBASE'] = self.gff.getEndBase( gene )
        pe['INTRON'] = self.gff.get_non_cds( gene )                
        pe['PRODUCT-TYPE'] = self.annot.getProductType( geneId )
        name = self.annot.getName( geneId )
        if name == geneId:
            pe['NAME'] = pe['ID']
        else:
            pe['NAME'] = name
        pe['DBLINK']  = self.annot.getDBlinks( geneId )
        pe['functions'] = self.annot.getFunctions( geneId )
        return pe
    
    def get_entries( self, geneKey ):
        pathologic_entries = {}
        for gene in self.gff.getGenes():
            pe = self.gene2entry( gene, geneKey )
            if pe['genetic_element'] in pathologic_entries:
                pathologic_entries[pe['genetic_element']].append( pe )
            else:
                pathologic_entries[pe['genetic_element']] = [pe]
        return pathologic_entries
    def generate_pathologic_files( self, pathologic_entries, pathologic_files={},template_file='genetic_element_{}.pf' , outputdir='.'):
        for genetic_element in pathologic_entries:
            if genetic_element not in pathologic_files:
                pathologic_files[genetic_element] = template_file.format(genetic_element )
                with open(os.path.join(outputdir,pathologic_files[genetic_element]), 'w') as pf:
                    for pe in sorted(pathologic_entries[genetic_element], key=operator.itemgetter('ID')):
                        pf.write(str(PathologicEntry(**pe)))
            else:
                with open(os.path.join(outputdir,pathologic_files[pe['genetic_element']]),'a') as pf:
                    for pe in sorted(pathologic_entries[genetic_element], key=operator.itemgetter('ID')):
                        pf.write(str(PathologicEntry(**pe)))
        return pathologic_files

    def generate_genetic_elements_file( self, pathologic_files,seq_files={}, circular={},outputdir='.'):
        with open(os.path.join(outputdir,'genetic-elements.dat'),'w') as out:
            for genetic_element in sorted(pathologic_files):
                out.write('ID\t{}\n'.format(genetic_element))
                out.write('NAME\t{}\n'.format( genetic_element))
                out.write('TYPE\t{}\n'.format( ':CHRSM'))
                if genetic_element in circular:
                    out.write('CIRCULAR?\t{}\n'.format(circular[genetic_element]))
                else:
                    out.write('CIRCULAR\tN\n')
                out.write('ANNOT-FILE\t{}\n'.format(pathologic_files[genetic_element]))
                if genetic_element in seq_files:
                    out.write('SEQ-FILE\t{}\n'.format(seq_files[genetic_element]))
                out.write('//\n')

def df_to_dol( df, keycol, valuecol, newcol, sep='|'):
    dol = {}
    df.columns = [c.strip('#') for c in df.columns]
    for i in df.index:
        key = df.loc[i,keycol]
        value = df.loc[i, valuecol]
        if key in dol:
            dol[key] += '{}{}'.format(sep,value)
        else:
            dol[key] = str(value)
    return pd.Series(dol).to_frame(newcol)

def writeable_dir(prospective_dir):
  if not os.path.isdir(prospective_dir):
       os.mkdir(prospective_dir)
  if os.access(prospective_dir, os.R_OK) and os.access(prospective_dir, os.W_OK):
      return prospective_dir
  elif  os.access(prospective_dir, os.R_OK):
      raise argparse.ArgumentTypeError("{0} is not a writeable dir".format(prospective_dir))
  elif os.access(prospective_dir, os.W_OK):
      raise argparse.ArgumentTypeError("{0} is not a readable dir".format(prospective_dir))
  else:
      raise argparse.ArgumentTypeError("{0} is neither a readable nor a writeable dir".format(prospective_dir))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert annotations to Pathologic file format')
    parser.add_argument('--gff' , type=argparse.FileType('r'), help='input gff file name')
    parser.add_argument('--gtf' , type=argparse.FileType('r'), help='input gtf file name')
    parser.add_argument('--ec', type=argparse.FileType('r'), help='input EC number annotations file name')
    parser.add_argument('--kog', type=argparse.FileType('r'), help='input kog annotations file name')
    parser.add_argument('--go',type=argparse.FileType('r'), help='input go annotations file name')
    parser.add_argument('--seq', type=argparse.FileType('r'), help='input (unmasked) sequence file')
    parser.add_argument('--outputdir', type=writeable_dir, help='output directory')
    args = parser.parse_args()
    go = df_to_dol(pd.read_table(args.go.name), 'proteinId', 'goAcc', 'GO')
    godef = df_to_dol(pd.read_table(args.go.name), 'proteinId', 'goName', 'GOdef')
    kog = df_to_dol(pd.read_table(args.kog.name), 'proteinId', 'kogid', 'kog')
    kogdef = df_to_dol(pd.read_table(args.kog.name), 'proteinId', 'kogdefline', 'kogdef')
    ec = df_to_dol(pd.read_table(args.ec.name), 'proteinId', 'ecNum', 'EC')
    ecdef = df_to_dol(pd.read_table(args.ec.name), 'proteinId', 'definition', 'ECdef')

    annot = pd.concat([ec, go, kog, ecdef, godef, kogdef], axis=1, join='outer')
    if args.gff:
        g2p = GFFandAnnot2Pathologic(GFF2Pathologic(args.gff.name),JGIAnnot( annot  ))
    elif args.gtf:
        g2p = GFFandAnnot2Pathologic(GTF2Pathologic(args.gff.name),JGIAnnot( annot  ))
    pe = g2p.get_entries('proteinId')
    pf_files = g2p.generate_pathologic_files(pe, {}, '{}.pf', args.outputdir)
    seq_files = dict([(ge,'{}.fna'.format(ge)) for ge in pf_files])
    g2p.generate_genetic_elements_file(pf_files, seq_files, {}, args.outputdir)
    with open(args.seq.name,'r') as assembly:
        for record in SeqIO.parse(assembly, "fasta"):
            with open(os.path.join(args.outputdir,"{}.fna".format(record.id)), 'w') as out:
                SeqIO.write([record], out, "fasta")
