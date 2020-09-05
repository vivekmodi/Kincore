#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys
from flask import Flask, render_template, url_for, redirect, request
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from sqlalchemy import func
from flask import Markup
#from flask_wtf import FlaskForm
#from flask import session
#from wtforms import StringField, SubmitField
#from wtforms.validators import DataRequired
#from werkzeug.utils import secure_filename
#from collections import defaultdict
from sqlalchemy import asc


pwd=os.getcwd()
#pwd='/var/www/html/site'
sys.path.append(pwd+'/scripts')        #To import webserver_script from this directory
from webserver_script import identify_state

UPLOAD_FOLDER = (pwd+'/server/uploads')
ALLOWED_EXTENSIONS = {'gz', 'pdb','cif'}

app=Flask(__name__)
app.config['SECRET_KEY'] = 'mysecretkey'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['ALLOWED_EXTENSIONS'] = ALLOWED_EXTENSIONS
app.jinja_env.add_extension('jinja2.ext.loopcontrols')

############SQL Databse and Models############################
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///data_new.sqlite'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

db = SQLAlchemy(app)
Migrate(app,db)

class Cluster(db.Model):
    pdb=db.Column(db.Text,primary_key=True)
    protein_name=db.Column(db.Text)
    uniprotid=db.Column(db.Text)
    chain_length=db.Column(db.Integer)
    str_begin=db.Column(db.Integer)
    str_end=db.Column(db.Integer)
    specie=db.Column(db.Text)
    resolution=db.Column(db.REAL)
    method=db.Column(db.Text)
    rvalue=db.Column(db.REAL)
    freervalue=db.Column(db.REAL)
    uniseq=db.Column(db.Text)
    gene=db.Column(db.Text)
    domain=db.Column(db.Text)
    group=db.Column(db.Text)
    uniprotacc=db.Column(db.Text)
    domain_begin=db.Column(db.Integer)
    domain_end=db.Column(db.Integer)
    alkres=db.Column(db.Text)
    alknum=db.Column(db.Integer)
    rreres=db.Column(db.Text)
    rrenum=db.Column(db.Integer)
    hrdres=db.Column(db.Text)
    hrdnum=db.Column(db.Integer)
    dfgres=db.Column(db.Text)
    dfgnum=db.Column(db.Integer)
    aperes=db.Column(db.Text)
    apenum=db.Column(db.Integer)
    dfgaspres=db.Column(db.Text)
    gtknum=db.Column(db.Integer)
    gtkres=db.Column(db.Text)
    hinge1=db.Column(db.Text)
    status=db.Column(db.Text)
    x_o_edia=db.Column(db.REAL)
    asp_o_edia=db.Column(db.REAL)
    phe_o_edia=db.Column(db.REAL)
    gly_o_edia=db.Column(db.REAL)
    x_phi=db.Column(db.REAL)
    x_psi=db.Column(db.REAL)
    asp_phi=db.Column(db.REAL)
    asp_psi=db.Column(db.REAL)
    asp_chi1=db.Column(db.REAL)
    asp_chi2=db.Column(db.REAL)
    phe_phi=db.Column(db.REAL)
    phe_psi=db.Column(db.REAL)
    phe_chi1=db.Column(db.REAL)
    phe_chi2=db.Column(db.REAL)
    gly_phi=db.Column(db.REAL)
    gly_psi=db.Column(db.REAL)
    domainBreak=db.Column(db.Integer)
    loopBreak=db.Column(db.Integer)
    synonym=db.Column(db.Text)
    x_mut=db.Column(db.Text)
    asp_mut=db.Column(db.Text)
    phe_mut=db.Column(db.Text)
    gly_mut=db.Column(db.Text)
    chain_mut=db.Column(db.Text)
    chain_phos=db.Column(db.Text)
    modified_aa=db.Column(db.Text)
    first_obs_res=db.Column(db.Integer)
    last_obs_res=db.Column(db.Integer)
    xdfg_resolved=db.Column(db.Text)
    ligand=db.Column(db.Text)
    strseq=db.Column(db.Text)
    seqres=db.Column(db.Text)
    lys_glu=db.Column(db.REAL)
    phe_glu4=db.Column(db.REAL)
    phe_lys=db.Column(db.REAL)
    chelix=db.Column(db.Text)
    spatial=db.Column(db.Text)
    dihedral=db.Column(db.Text)
    color=db.Column(db.Text)
    ligand_type=db.Column(db.Text)


    def __init__(self,pdb,protein_name,uniprotid,chain_length,str_begin,str_end,specie,resolution,method,rvalue,freervalue,uniseq,gene,\
    domain,group,uniprotacc,domain_begin,domain_end,alkres,alknum,rreres,rrenum,hrdres,hrdnum,dfgres,dfgnum,aperes,apenum,dfgaspres,\
    gtknum,gtkres,hinge1,status,x_o_edia,asp_o_edia,phe_o_edia,gly_o_edia,x_phi,x_psi,asp_phi,asp_psi,\
    asp_chi1,asp_chi2,phe_phi,phe_psi,phe_chi1,phe_chi2,gly_phi,gly_psi,domainBreak,loopBreak,synonym,x_mut,asp_mut,phe_mut,gly_mut,chain_mut,\
    chain_phos,modified_aa,first_obs_res,last_obs_res,ligand,xdfg_resolved,strseq,seqres,lys_glu,phe_glu4,phe_lys,chelix,spatial,dihedral,\
    color,ligand_type):
        self.pdb=pdb;self.protein_name=protein_name;self.uniprotid=uniprotid;self.chain_length=chain_length;self.str_begin=str_begin;\
        self.str_end=str_end;self.specie=specie;self.resolution=resolution;self.method=method;self.rvalue=rvalue;self.freervalue=freervalue;\
        self.uniseq=uniseq;self.gene=gene;self.domain=domain;self.group=group;self.uniprotacc=uniprotacc;self.domain_begin=domain_begin;\
        self.domain_end=domain_end;self.alkres=alkres;self.alknum=alknum;self.rreres=rreres;self.rrenum=rrenum;self.hrdres=hrdres;self.hrdnum=hrdnum;\
        self.dfgres=dfgres;self.dfgnum=dfgnum;self.aperes=aperes;self.apenum=apenum;self.dfgaspres=dfgaspres;self.gtknum=gtknum;\
        self.gtkres=gtkres;self.hinge1=hinge1;self.status=status;self.x_o_edia=x_o_edia;\
        self.asp_o_edia=asp_o_edia;self.phe_o_edia=phe_o_edia;self.gly_o_edia=gly_o_edia;\
        self.x_phi=x_phi;self.x_psi=x_psi;self.asp_phi=asp_phi;self.asp_psi=asp_psi;self.asp_chi1=asp_chi1;self.asp_chi2=asp_chi2;self.phe_phi=phe_phi;\
        self.phe_psi=phe_psi;self.phe_chi1=phe_chi1;self.phe_chi2=phe_chi2;self.gly_phi=gly_phi;self.gly_psi=gly_psi;self.domainBreak=domainBreak;\
        self.loopBreak=loopBreak;self.synonym=synonym;self.x_mut=x_mut;self.asp_mut=asp_mut;self.phe_mut=phe_mut;self.gly_mut=gly_mut;\
        self.chain_mut=chain_mut;self.chain_phos=chain_phos;self.modified_aa=modified_aa;self.first_obs_res=first_obs_res;self.last_obs_res=last_obs_res;\
        self.ligand=ligand;self.xdfg_resolved=xdfg_resolved;self.strseq=strseq;self.seqres=seqres;self.lys_glu=lys_glu;\
        self.phe_glu4=phe_glu4;self.phe_lys=phe_lys;self.chelix=chelix;self.spatial=spatial;self.dihedral=dihedral;\
        self.color=color;self.ligand_type=ligand_type;


    def __repr__(self):
        return f'{self.pdb} {self.protein_name} {self.uniprotid} {self.chain_length} {self.str_begin} {self.str_end} {self.specie} {self.resolution}\
        {self.method} {self.rvalue} {self.freervalue} {self.uniseq} {self.gene} {self.domain} {self.group} {self.uniprotacc} {self.domain_begin} {self.domain_end}\
        {self.alkres} {self.alknum} {self.rreres} {self.rrenum} {self.hrdres} {self.hrdnum} {self.dfgres} {self.dfgnum} {self.aperes}\
        {self.apenum} {self.dfgaspres} {self.gtknum} {self.gtkres} {self.hinge1} {self.status} {self.x_o_edia} \
        {self.asp_o_edia} {self.phe_o_edia} {self.gly_o_edia} {self.x_phi} {self.x_psi} {self.asp_phi} {self.asp_psi}\
        {self.asp_chi1} {self.asp_chi2} {self.phe_phi} {self.phe_psi} {self.phe_chi1} {self.phe_chi2} {self.gly_phi} {self.gly_psi} {self.domainBreak}\
        {self.loopBreak} {self.synonym} {self.x_mut} {self.asp_mut} {self.phe_mut} {self.gly_mut} {self.chain_mut} {self.chain_phos} {self.modified_aa}\
        {self.first_obs_res} {self.last_obs_res} {self.ligand} {self.xdfg_resolved} {self.strseq} {self.seqres} {self.lys_glu}\
        {self.phe_glu4} {self.phe_lys} {self.chelix} {self.spatial} {self.dihedral} {self.color} {self.ligand_type}'


###################FUNCTIONS###############################
def create_lists():
    pdbListDb=list();chainListDb=list();uniprotIdListDb=list();uniprotAccListDb=list();geneListDb=list();ligandListDb=list();specieListDb=list()
    for pdbs in Cluster.query.with_entities(Cluster.pdb):
        pdbListDb.append(str(pdbs[0])[0:4])
    for pdbs in Cluster.query.with_entities(Cluster.pdb):
        chainListDb.append(pdbs[0])
    for uniprots in Cluster.query.with_entities(Cluster.uniprotid):
        uniprotIdListDb.append(uniprots[0])
    for uniprots in Cluster.query.with_entities(Cluster.uniprotacc):
        uniprotAccListDb.append(uniprots[0])
    for genes in Cluster.query.with_entities(Cluster.gene):
        geneListDb.append(genes[0])
    for ligands in Cluster.query.with_entities(Cluster.ligand):
        ligandListDb.append(ligands[0])
    for species in Cluster.query.with_entities(Cluster.specie):
        specieListDb.append(species[0])

    return (pdbListDb,chainListDb,uniprotIdListDb,uniprotAccListDb,geneListDb,ligandListDb,specieListDb)

def create_color_lists():
    clusterColor={'BLAminus':'rgb(179,215,229)','BLAplus':'rgb(254, 214, 154)','ABAminus':'rgba(243,140,184)','BLBminus':'rgb(250,128,114)',\
    'BLBplus':'rgb(199,233,175)','BLBtrans':'rgb(253,156,104)','BABtrans':'rgb(2,71,254)','BBAminus':'rgb(154,180,254)','None':'rgb(230,230,230)'}
    return clusterColor

def count_structures_groups(groupList):       #given a list this functions counts the number of chains in each conformation
    strCount=dict();geneList=list();totalGroup=dict();totalDihedral=dict()
    groups=('AGC','CAMK','CK1','CMGC','NEK','RGC','STE','TKL','TYR','OTHER')
    spatialList=('DFGin','DFGinter','DFGout','None')
    dihedralList=('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','BABtrans','BBAminus','None')

    for group in groups:
        for spatial in spatialList:
            for dihedral in dihedralList:
                strCount[group,spatial,dihedral]=0
                totalDihedral[spatial,dihedral]=0
                totalGroup[group]=0

    for item in groupList:
        geneList.append(item.gene)
        strCount[item.group,item.spatial,item.dihedral]+=1
        totalDihedral[item.spatial,item.dihedral]+=1
        totalGroup[item.group]+=1
    geneCount=len(set(geneList))
    return (strCount,geneCount,totalGroup,totalDihedral)

def count_structures_all(groupList):       #given a list this function counts the number of chains in each conformation
    strCount=dict();geneList=list()
    spatialList=('DFGin','DFGinter','DFGout','None')
    dihedralList=('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','BABtrans','BBAminus','None')

    for spatial in spatialList:
        for dihedral in dihedralList:
            strCount[spatial,dihedral]=0

    for item in groupList:
        geneList.append(item.gene)
        strCount[item.spatial,item.dihedral]+=1
    geneCount=len(set(geneList))
    return (strCount,geneCount)

def min_atom_missing(group_list):
    subListRepr=dict()
    for item in group_list:
        subListRepr[item.pdb]=(int(item.domainBreak),float(item.resolution))

    count=0
    for entry in sorted(subListRepr,key=subListRepr.get):
        if count==0:
            count=1
            return entry              #just return the first entry from sorted list
    return 'None'

def write_text_file(sublist,csvFile):
    fhandle_textFile=open(f'{pwd}/static/{csvFile}','w')
    fhandle_textFile.write('Organism\tGroup\tGene\tUniprotID\tPDB\tMethod\tResolution\tRfac\tFreeRfac\tSpatialLabel\tDihedralLabel\tC-helix\tLigand\tLigandType\tDFG_Phe\tEdia_X_O\tEdia_Asp_O\tEdia_Phe_O\tEdia_Gly_O\tProteinName\n')

    for item in sublist:
        ligand_list=list()
        if ',' in item.ligand:
            lig_items=item.ligand.split(',')
            for lig in lig_items:
                ligandname=lig.split(':')[0];ligand_list.append(ligandname)
        else:
            ligandname=item.ligand.split(':')[0]
            ligand_list.append(ligandname)
        ligand_list=','.join(ligand_list)
        fhandle_textFile.write(f'{item.specie}\t{item.group}\t{item.domain}\t{item.uniprotid}\t{item.pdb}\t{item.method}\t{item.resolution}\t{item.rvalue}\t{item.freervalue}\t{item.spatial}\t{item.dihedral}\t{item.chelix}\t{ligand_list}\t{item.ligand_type}\t{item.dfgnum}\t{item.x_o_edia}\t{item.asp_o_edia}\t{item.phe_o_edia}\t{item.gly_o_edia}\t{item.protein_name}\n')
    fhandle_textFile.close()

##########ROUTES##########################

@app.route('/index')
@app.route('/home')
@app.route('/')
def home():
    return render_template('home.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/browse')
def browse():
    return redirect(url_for('multipleQuery',groupSelect='All',labelSelect='All',ligTypeSelect='All'))

@app.route('/geneListHelp')
def geneListHelp():
	fhandle_genehelp=open('{pwd}/static/geneListHelpFile.tab','r')
	genehelp=fhandle_genehelp.read()
	fhandle_genehelp.close()
	return render_template('genelisthelp.html', text=genehelp)

@app.route('/formSearch', methods=['GET','POST'])
def formSearch():
    (pdbListDb,chainListDb,uniprotIdListDb,uniprotAccListDb,geneListDb,ligandListDb,specieListDb)=create_lists()

    if request.method=='POST':
        inputString=request.form['inputString'].upper()

        if inputString in pdbListDb:    #match without chain
            return redirect(url_for('uniqueQuery',queryname=inputString,settings='PDB'))
        if inputString in chainListDb:  #match with chain
            return redirect(url_for('uniqueQuery',queryname=inputString,settings='PDB'))
        if inputString in uniprotIdListDb:
            return redirect(url_for('uniqueQuery',queryname=inputString,settings='UNIPROTID'))
        if inputString in geneListDb:                 
            return redirect(url_for('uniqueQuery',queryname=inputString,settings='GENE'))
        for items in ligandListDb:          #Loop is required because some entries have two ligands
            if inputString in items:
                return redirect(url_for('uniqueQuery',queryname=inputString,settings='LIGAND'))
        else:
            return render_template('nomatch.html')

    return render_template ('search.html')

@app.route('/formSearchMultiple',methods=['GET','POST'])
def formSearchMultiple():
    if request.method=='POST':
        groupSelect=request.form['groupSelect']
        labelSelect=request.form['labelSelect']
        ligTypeSelect=request.form['ligTypeSelect']
        return redirect(url_for('multipleQuery',groupSelect=groupSelect,labelSelect=labelSelect,ligTypeSelect=ligTypeSelect))
    return render_template ('search.html')



def allowed_filename(userfilename):
    if not '.' in userfilename:         #To check if there is a dot in the filename
        return False
    ext=userfilename.split('.')[1]

    if ext.lower() in app.config['ALLOWED_EXTENSIONS']:
        return True
    else:
        return False



@app.route('/alignment')
def alignment():
    return render_template('alignment.html')

@app.route('/biojs')
def biojs():
    return render_template('Human-PK-BioJS.html')

@app.route('/phylogeny')
def phylogeny():
    return render_template('phylogeny.html')

@app.route('/download')
def download():
    pdb_list=Cluster.query.all()
    agcDomains=list();camkDomains=list();ck1Domains=list();cmgcDomains=list();nekDomains=list();steDomains=list();tklDomains=list();\
    tyrDomains=list();otherDomains=list();ligandList=list()
    for items in pdb_list:
        if items.specie!='Homo sapiens':
            continue
        if items.group=='AGC':
            agcDomains.append(items.domain)
        if items.group=='CAMK':
            camkDomains.append(items.domain)
        if items.group=='CK1':
            ck1Domains.append(items.domain)
        if items.group=='CMGC':
            cmgcDomains.append(items.domain)
        if items.group=='NEK':
            nekDomains.append(items.domain)
        if items.group=='STE':
            steDomains.append(items.domain)
        if items.group=='TKL':
            tklDomains.append(items.domain)
        if items.group=='TYR':
            tyrDomains.append(items.domain)
        if items.group=='OTHER':
            otherDomains.append(items.domain)
        if ',' in items.ligand:
            lig_items=items.ligand.split(',')
            for lig in lig_items:
                ligandname=lig.split(':')[0]
                ligandList.append(ligandname)
        else:
            ligandname=items.ligand.split(':')[0]
            ligandList.append(ligandname)

    agcDomains=sorted(list(set(agcDomains)));camkDomains=sorted(list(set(camkDomains)));ck1Domains=sorted(list(set(ck1Domains)));\
    cmgcDomains=sorted(list(set(cmgcDomains)));nekDomains=sorted(list(set(nekDomains)));steDomains=sorted(list(set(steDomains)));\
    tklDomains=sorted(list(set(tklDomains)));tyrDomains=sorted(list(set(tyrDomains)));otherDomains=sorted(list(set(otherDomains)));
    ligandList=sorted(list(set(ligandList)))

    return render_template('download.html',agcDomains=agcDomains,camkDomains=camkDomains,ck1Domains=ck1Domains,cmgcDomains=cmgcDomains,\
    nekDomains=nekDomains,steDomains=steDomains,tklDomains=tklDomains,tyrDomains=tyrDomains,otherDomains=otherDomains,ligandList=ligandList)

@app.route('/contact')
def contact():
    return render_template('contact.html')


@app.route('/help')
def help():
    return render_template('help.html')



@app.route('/multipleQuery/<groupSelect>/<labelSelect>/<ligTypeSelect>')
def multipleQuery(groupSelect,labelSelect,ligTypeSelect):

    group_list=dict();total_count=dict();strCount=dict();geneCount=dict();totalGroup=dict();totalDihedral=dict()
    gene_list=dict();total_count=dict();reprStr=dict();pymolSession=dict();pymolSessionRe=dict();pymolScript=dict();pymolScriptRe=dict();coordinateFiles=dict()
    nglList=dict();dfgNumReprStr=dict();subList=dict();subListPymol=dict();dfgNumReprStr['Human']=list();dfgNumReprStr['All']=list();dfgNumReprStr['Nonhuman']=list()
    csvFile=dict()
    clusterColor=create_color_lists()

    ligand_type=ligTypeSelect
    if ligTypeSelect=='All':
        ligTypeSelect=''
        ligand_type='All'
    if ligTypeSelect=='Type1':      #without this condition Type1 also matches Type1.5
        dontmatch='%Type1.5%'
    else:
        dontmatch=''
    if ligTypeSelect=='Type1.5_Front' or ligTypeSelect=='Type1.5_Back':
        ligTypeSelect='Type1.5';ligand_type='Type1.5'
    if groupSelect=='All':
        groupSelect=''

    if labelSelect=='All' and groupSelect=='':    #condition for All-All
        for tabs in ('Human','All','Nonhuman'):
            if tabs=='All':
                group_list[tabs]=Cluster.query.filter(Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
                total_count[tabs]=Cluster.query.filter(Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
                (strCount_all,geneCount_all,totalGroup_all,totalDihedral_all)=count_structures_groups(group_list[tabs])
            if tabs=='Human':
                group_list[tabs]=Cluster.query.filter(Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
                total_count[tabs]=Cluster.query.filter(Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
                (strCount_human,geneCount_human,totalGroup_human,totalDihedral_human)=count_structures_groups(group_list[tabs])
            if tabs=='Nonhuman':
                group_list[tabs]=Cluster.query.filter(Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
                total_count[tabs]=Cluster.query.filter(Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
                (strCount_nonhuman,geneCount_nonhuman,totalGroup_nonhuman,totalDihedral_nonhuman)=count_structures_groups(group_list[tabs])

            csvFile[tabs]=f'downloads/text-files/{tabs}_Allgroups_Allspatials_Alldihedrals_{ligand_type}.tab'
            write_text_file(group_list[tabs],csvFile[tabs])

        return render_template('All-All.html',group_list=group_list,total_count=total_count,strCount_all=strCount_all,strCount_human=strCount_human,strCount_nonhuman=strCount_nonhuman,\
        geneCount_all=geneCount_all,geneCount_human=geneCount_human,geneCount_nonhuman=geneCount_nonhuman,totalGroup_all=totalGroup_all,totalGroup_human=totalGroup_human,\
        totalGroup_nonhuman=totalGroup_nonhuman,totalDihedral_all=totalDihedral_all,totalDihedral_human=totalDihedral_human,totalDihedral_nonhuman=totalDihedral_nonhuman,\
        csvFile=csvFile,ligand_type=ligand_type)

    if  groupSelect=='' and labelSelect in ('DFGin','DFGinter','DFGout','None'):    #condition for All-spatial

        for tabs in ('Human','All','Nonhuman'):
            group_list['All']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['All']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_all,geneCount_all,totalGroup_all,totalDihedral_all)=count_structures_groups(group_list['All'])

            group_list['Human']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Human']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_human,geneCount_human,totalGroup_human,totalDihedral_human)=count_structures_groups(group_list['Human'])

            group_list['Nonhuman']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Nonhuman']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_nonhuman,geneCount_nonhuman,totalGroup_nonhuman,totalDihedral_nonhuman)=count_structures_groups(group_list['Nonhuman'])

            if labelSelect=='DFGin':
                 for dihedral in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','None'):
                     subList['Human']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['All']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['Nonhuman']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()

                     reprStr[tabs,labelSelect,dihedral]=min_atom_missing(subList[tabs])
                     nglList[tabs,reprStr[tabs,labelSelect,dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).all()

                     if Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).first():
                         dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).first().dfgnum)

            if labelSelect=='DFGinter':
                 for dihedral in ('BABtrans','None'):
                     subList['Human']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['All']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['Nonhuman']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()

                     reprStr[tabs,labelSelect,dihedral]=min_atom_missing(subList[tabs])
                     nglList[tabs,reprStr[tabs,labelSelect,dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).all()

                     if Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).first():
                         dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).first().dfgnum)

            if labelSelect=='DFGout':
                 for dihedral in ('BBAminus','None'):
                     subList['Human']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['All']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['Nonhuman']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()

                     reprStr[tabs,labelSelect,dihedral]=min_atom_missing(subList[tabs])
                     nglList[tabs,reprStr[tabs,labelSelect,dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).all()
                     if Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).first():
                         dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).first().dfgnum)

            if labelSelect=='None':
                     subList['Human']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['All']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['Nonhuman']=Cluster.query.filter(Cluster.spatial==labelSelect,Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()

                     reprStr[tabs,'None','None']=min_atom_missing(subList[tabs])
                     nglList[tabs,reprStr[tabs,'None','None']]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).all()
                     if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).first():
                         dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).first().dfgnum)

            pymolSession[tabs]=f'downloads/pymolSessions/{tabs}_Allgroups_{labelSelect}_Alldihedrals_{ligand_type}.pse.zip'
            pymolScript[tabs]=f'downloads/pymolSessionScripts/{tabs}_Allgroups_{labelSelect}_Alldihedrals_{ligand_type}.zip'
            pymolSessionRe[tabs]=f'downloads/pymolSessions/Repr_{tabs}_Allgroups_{labelSelect}_Alldihedrals_{ligand_type}.pse.zip'
            pymolScriptRe[tabs]=f'downloads/pymolSessionScripts/Repr_{tabs}_Allgroups_{labelSelect}_Alldihedrals_{ligand_type}.zip'
            coordinateFiles[tabs]=f'downloads/coordinateFiles/{tabs}_Allgroups_{labelSelect}_Alldihedrals_{ligand_type}'

            csvFile[tabs]=f'downloads/text-files/{tabs}_Allgroups_{labelSelect}_Alldihedrals_{ligand_type}.tab'
            write_text_file(group_list[tabs],csvFile[tabs])

        return render_template('All-Spatial.html',group_list=group_list,total_count=total_count,strCount_human=strCount_human,geneCount_human=geneCount_human,\
        strCount_all=strCount_all,geneCount_all=geneCount_all,strCount_nonhuman=strCount_nonhuman,\
        geneCount_nonhuman=geneCount_nonhuman,reprStr=reprStr,label=labelSelect,pymolSession=pymolSession,pymolSessionRe=pymolSessionRe,pymolScript=pymolScript,\
        pymolScriptRe=pymolScriptRe,clusterColor=clusterColor,coordinateFiles=coordinateFiles,csvFile=csvFile,\
        nglList=nglList,dfgNumReprStr=dfgNumReprStr,ligand_type=ligand_type)

    if  groupSelect=='' and labelSelect in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','BABtrans','BBAminus'):   #condition for All-dihedral

        for tabs in ('Human','All','Nonhuman'):
            group_list['All']=Cluster.query.filter(Cluster.dihedral==labelSelect,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['All']=Cluster.query.filter(Cluster.dihedral==labelSelect,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_all,geneCount_all,totalGroup_all,totalDihedral_all)=count_structures_groups(group_list['All'])

            group_list['Human']=Cluster.query.filter(Cluster.dihedral==labelSelect,Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Human']=Cluster.query.filter(Cluster.dihedral==labelSelect,Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_human,geneCount_human,totalGroup_human,totalDihedral_human)=count_structures_groups(group_list['Human'])

            group_list['Nonhuman']=Cluster.query.filter(Cluster.dihedral==labelSelect,Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Nonhuman']=Cluster.query.filter(Cluster.dihedral==labelSelect,Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_nonhuman,geneCount_nonhuman,totalGroup_nonhuman,totalDihedral_nonhuman)=count_structures_groups(group_list['Nonhuman'])

            dihedral=labelSelect
            subList['Human']=Cluster.query.filter(Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
            subList['All']=Cluster.query.filter(Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
            subList['Nonhuman']=Cluster.query.filter(Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()

            reprStr[tabs,'All',dihedral]=min_atom_missing(subList[tabs])
            nglList[tabs,reprStr[tabs,'All',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'All',dihedral],Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
            if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'All',dihedral]).first():
                dfgNumReprStr[tabs]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'All',dihedral],Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).first().dfgnum

            if labelSelect in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans'):
                spatial='DFGin'
                pymolSession[tabs]=f'downloads/pymolSessions/{tabs}_Allgroups_DFGin_{labelSelect}_{ligand_type}.pse.zip'
                pymolScript[tabs]=f'downloads/pymolSessionScripts/{tabs}_Allgroups_DFGin_{labelSelect}_{ligand_type}.zip'
                pymolSessionRe[tabs]=f'downloads/pymolSessions/Repr_{tabs}_Allgroups_DFGin_{labelSelect}_{ligand_type}.pse.zip'
                pymolScriptRe[tabs]=f'downloads/pymolSessionScripts/Repr_{tabs}_Allgroups_DFGin_{labelSelect}_{ligand_type}.zip'
                csvFile[tabs]=f'downloads/text-files/{tabs}_Allgroups_DFGin_{labelSelect}_{ligand_type}.tab'
                write_text_file(group_list[tabs],csvFile[tabs])
                coordinateFiles[tabs]=f'downloads/coordinateFiles/{tabs}_Allgroups_DFGin_{labelSelect}_{ligand_type}'
            if labelSelect=='BABtrans':
                spatial='DFGinter'
                pymolSession[tabs]=f'downloads/pymolSessions/{tabs}_Allgroups_DFGinter_{labelSelect}_{ligand_type}.pse.zip'
                pymolScript[tabs]=f'downloads/pymolSessionScripts/{tabs}_Allgroups_DFGinter_{labelSelect}_{ligand_type}.zip'
                pymolSessionRe[tabs]=f'downloads/pymolSessions/Repr_{tabs}_Allgroups_DFGinter_{labelSelect}_{ligand_type}.pse.zip'
                pymolScriptRe[tabs]=f'downloads/pymolSessionScripts/Repr_{tabs}_Allgroups_DFGinter_{labelSelect}_{ligand_type}.zip'
                csvFile[tabs]=f'downloads/text-files/{tabs}_Allgroups_DFGinter_{labelSelect}_{ligand_type}.tab'
                write_text_file(group_list[tabs],csvFile[tabs])
                coordinateFiles[tabs]=f'downloads/coordinateFiles/{tabs}_Allgroups_DFGinter_{labelSelect}_{ligand_type}'
            if labelSelect=='BBAminus':
                spatial='DFGout'
                pymolSession[tabs]=f'downloads/pymolSessions/{tabs}_Allgroups_DFGout_{labelSelect}_{ligand_type}.pse.zip'
                pymolScript[tabs]=f'downloads/pymolSessionScripts/{tabs}_Allgroups_DFGout_{labelSelect}_{ligand_type}.zip'
                pymolSessionRe[tabs]=f'downloads/pymolSessions/Repr_{tabs}_Allgroups_DFGout_{labelSelect}_{ligand_type}.pse.zip'
                pymolScriptRe[tabs]=f'downloads/pymolSessionScripts/Repr_{tabs}_Allgroups_DFGout_{labelSelect}_{ligand_type}.zip'
                csvFile[tabs]=f'downloads/text-files/{tabs}_Allgroups_DFGout_{labelSelect}_{ligand_type}.tab'
                write_text_file(group_list[tabs],csvFile[tabs])
                coordinateFiles[tabs]=f'downloads/coordinateFiles/{tabs}_Allgroups_DFGout_{labelSelect}_{ligand_type}'

        return render_template('All-dihedral.html',group_list=group_list,total_count=total_count,spatial=spatial,strCount_human=strCount_human,strCount_all=strCount_all,\
        strCount_nonhuman=strCount_nonhuman,reprStr=reprStr,geneCount_human=geneCount_human,geneCount_all=geneCount_all,geneCount_nonhuman=geneCount_nonhuman,\
        label=labelSelect,pymolSession=pymolSession,pymolScript=pymolScript,pymolSessionRe=pymolSessionRe,pymolScriptRe=pymolScriptRe,clusterColor=clusterColor,\
        coordinateFiles=coordinateFiles,csvFile=csvFile,nglList=nglList,dfgNumReprStr=dfgNumReprStr,ligand_type=ligand_type)

    if groupSelect!='' and labelSelect=='All':   #condition for Group-All

        for tabs in ('Human','All','Nonhuman'):
            group_list['All']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['All']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_all,geneCount_all,totalGroup_all,totalDihedral_all)=count_structures_groups(group_list['All'])

            group_list['Human']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Human']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_human,geneCount_human,totalGroup_human,totalDihedral_human)=count_structures_groups(group_list['Human'])

            group_list['Nonhuman']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Nonhuman']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_nonhuman,geneCount_nonhuman,totalGroup_nonhuman,totalDihedral_nonhuman)=count_structures_groups(group_list['Nonhuman'])

            for dihedral in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','None'):
                subList['Human']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie=='Homo sapiens',Cluster.spatial=='DFGin',Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                subList['All']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.spatial=='DFGin',Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                subList['Nonhuman']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie!='Homo sapiens',Cluster.spatial=='DFGin',Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                reprStr[tabs,'DFGin',dihedral]=min_atom_missing(subList[tabs])
                nglList[tabs,reprStr[tabs,'DFGin',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGin',dihedral]).all()
                if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGin',dihedral]).first():
                    dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGin',dihedral]).first().dfgnum)

            for dihedral in ('BABtrans','None'):
                subList['Human']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie=='Homo sapiens',Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                subList['All']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                subList['Nonhuman']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie!='Homo sapiens',Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                reprStr[tabs,'DFGinter',dihedral]=min_atom_missing(subList[tabs])
                nglList[tabs,reprStr[tabs,'DFGinter',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGinter',dihedral]).all()
                if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGinter',dihedral]).first():
                    dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGinter',dihedral]).first().dfgnum)

            for dihedral in ('BBAminus','None'):
                subList['Human']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie=='Homo sapiens',Cluster.spatial=='DFGout',Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                subList['All']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.spatial=='DFGout',Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                subList['Nonhuman']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie!='Homo sapiens',Cluster.spatial=='DFGout',Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                reprStr[tabs,'DFGout',dihedral]=min_atom_missing(subList[tabs])
                nglList[tabs,reprStr[tabs,'DFGout',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGout',dihedral]).all()
                if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGout',dihedral]).first():
                    dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGout',dihedral]).first().dfgnum)

            subList['Human']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie=='Homo sapiens',Cluster.spatial=='None',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
            subList['All']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.spatial=='None',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
            subList['Nonhuman']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie!='Homo sapiens',Cluster.spatial=='None',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
            reprStr[tabs,'None','None']=min_atom_missing(subList[tabs])
            nglList[tabs,reprStr[tabs,'None','None']]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).all()
            if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).first():
                dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).first().dfgnum)

            pymolSession[tabs]=f'downloads/pymolSessions/{tabs}_{groupSelect}_Allspatials_Alldihedrals_{ligand_type}.pse.zip'
            pymolScript[tabs]=f'downloads/pymolSessionScripts/{tabs}_{groupSelect}_Allspatials_Alldihedrals_{ligand_type}.zip'
            pymolSessionRe[tabs]=f'downloads/pymolSessions/Repr_{tabs}_{groupSelect}_Allspatials_Alldihedrals_{ligand_type}.pse.zip'
            pymolScriptRe[tabs]=f'downloads/pymolSessionScripts/Repr_{tabs}_{groupSelect}_Allspatials_Alldihedrals_{ligand_type}.zip'

            csvFile[tabs]=f'downloads/text-files/{tabs}_{groupSelect}_Allspatials_Alldihedrals_{ligand_type}.tab'
            write_text_file(group_list[tabs],csvFile[tabs])
            coordinateFiles[tabs]=f'downloads/coordinateFiles/{tabs}_{groupSelect}_Allspatials_Alldihedrals_{ligand_type}'

        return render_template('Group-All.html',group_list=group_list,total_count=total_count,strCount_human=strCount_human,strCount_all=strCount_all,\
        strCount_nonhuman=strCount_nonhuman,reprStr=reprStr,geneCount_human=geneCount_human,geneCount_all=geneCount_all,geneCount_nonhuman=geneCount_nonhuman,\
        groupname=groupSelect,label=labelSelect,pymolSession=pymolSession,pymolScript=pymolScript,pymolSessionRe=pymolSessionRe,pymolScriptRe=pymolScriptRe,\
        clusterColor=clusterColor,coordinateFiles=coordinateFiles,csvFile=csvFile,nglList=nglList,dfgNumReprStr=dfgNumReprStr,ligand_type=ligand_type)

    if groupSelect!='' and labelSelect in ('DFGin','DFGinter','DFGout','None'):   #Group-Spatial

        for tabs in ('Human','All','Nonhuman'):
            group_list['All']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.spatial==labelSelect,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['All']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.spatial==labelSelect,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_all,geneCount_all,totalGroup_all,totalDihedral_all)=count_structures_groups(group_list['All'])

            group_list['Human']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.spatial==labelSelect,Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Human']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.spatial==labelSelect,Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_human,geneCount_human,totalGroup_human,totalDihedral_human)=count_structures_groups(group_list['Human'])

            group_list['Nonhuman']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.spatial==labelSelect,Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Nonhuman']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.spatial==labelSelect,Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_nonhuman,geneCount_nonhuman,totalGroup_nonhuman,totalDihedral_nonhuman)=count_structures_groups(group_list['Nonhuman'])

            if labelSelect=='DFGin':
                 for dihedral in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','None'):
                     subList['All']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['Human']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie=='Homo sapiens',Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['Nonhuman']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie!='Homo sapiens',Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()

                     reprStr[tabs,labelSelect,dihedral]=min_atom_missing(subList[tabs])
                     nglList[tabs,reprStr[tabs,labelSelect,dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).all()
                     if Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).first():
                         dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).first().dfgnum)

            if labelSelect=='DFGinter':
                 for dihedral in ('BABtrans','None'):
                     subList['All']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['Human']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie=='Homo sapiens',Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['Nonhuman']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie!='Homo sapiens',Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()

                     reprStr[tabs,labelSelect,dihedral]=min_atom_missing(subList[tabs])
                     nglList[tabs,reprStr[tabs,labelSelect,dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).all()
                     if Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).first():
                         dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).first().dfgnum)

            if labelSelect=='DFGout':
                 for dihedral in ('BBAminus','None'):
                     subList['All']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['Human']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie=='Homo sapiens',Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                     subList['Nonhuman']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie!='Homo sapiens',Cluster.spatial==labelSelect,Cluster.dihedral==dihedral,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()

                     reprStr[tabs,labelSelect,dihedral]=min_atom_missing(subList[tabs])
                     nglList[tabs,reprStr[tabs,labelSelect,dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).all()
                     if Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).first():
                         dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect,dihedral]).first().dfgnum)

            if labelSelect=='None':   #spatial None
                subList['All']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.spatial=='None',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                subList['Human']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie=='Homo sapiens',Cluster.spatial=='None',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
                subList['Nonhuman']=Cluster.query.filter(Cluster.group==groupSelect,Cluster.specie!='Homo sapiens',Cluster.spatial=='None',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()

                reprStr[tabs,'None','None']=min_atom_missing(subList[tabs])
                nglList[tabs,reprStr[tabs,'None','None']]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).all()
                if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).first():
                    dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).first().dfgnum)

            pymolSession[tabs]=f'downloads/pymolSessions/{tabs}_{groupSelect}_{labelSelect}_Alldihedrals_{ligand_type}.pse.zip'
            pymolScript[tabs]=f'downloads/pymolSessionScripts/{tabs}_{groupSelect}_{labelSelect}_Alldihedrals_{ligand_type}.zip'
            pymolSessionRe[tabs]=f'downloads/pymolSessions/Repr_{tabs}_{groupSelect}_{labelSelect}_Alldihedrals_{ligand_type}.pse.zip'
            pymolScriptRe[tabs]=f'downloads/pymolSessionScripts/Repr_{tabs}_{groupSelect}_{labelSelect}_Alldihedrals_{ligand_type}.zip'

            csvFile[tabs]=f'downloads/text-files/{tabs}_{groupSelect}_{labelSelect}_Alldihedrals_{ligand_type}.tab'
            write_text_file(group_list[tabs],csvFile[tabs])
            coordinateFiles[tabs]=f'downloads/coordinateFiles/{tabs}_{groupSelect}_{labelSelect}_Alldihedrals_{ligand_type}'

        return render_template('Group-Spatial.html',group_list=group_list,total_count=total_count,strCount_all=strCount_all,strCount_human=strCount_human,\
        strCount_nonhuman=strCount_nonhuman,reprStr=reprStr,geneCount_all=geneCount_all,geneCount_human=geneCount_human,geneCount_nonhuman=geneCount_nonhuman,\
        groupname=groupSelect,label=labelSelect,pymolSession=pymolSession,pymolScript=pymolScript,pymolSessionRe=pymolSessionRe,pymolScriptRe=pymolScriptRe,\
        clusterColor=clusterColor,coordinateFiles=coordinateFiles,csvFile=csvFile,nglList=nglList,dfgNumReprStr=dfgNumReprStr,ligand_type=ligand_type)

    if groupSelect!='' and labelSelect in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','BABtrans','BBAminus'):   #condition for Group-dihedral

        for tabs in ('Human','All','Nonhuman'):
            group_list['All']=Cluster.query.filter(Cluster.dihedral==labelSelect,Cluster.group.contains(groupSelect),Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['All']=Cluster.query.filter(Cluster.dihedral==labelSelect,Cluster.group.contains(groupSelect),Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_all,geneCount_all,totalGroup_all,totalDihedral_all)=count_structures_groups(group_list['All'])

            group_list['Human']=Cluster.query.filter(Cluster.dihedral==labelSelect,Cluster.group.contains(groupSelect),Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Human']=Cluster.query.filter(Cluster.dihedral==labelSelect,Cluster.group.contains(groupSelect),Cluster.specie=='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_human,geneCount_human,totalGroup_human,totalDihedral_human)=count_structures_groups(group_list['Human'])

            group_list['Nonhuman']=Cluster.query.filter(Cluster.dihedral==labelSelect,Cluster.group.contains(groupSelect),Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Nonhuman']=Cluster.query.filter(Cluster.dihedral==labelSelect,Cluster.group.contains(groupSelect),Cluster.specie!='Homo sapiens',Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()
            (strCount_nonhuman,geneCount_nonhuman,totalGroup_nonhuman,totalDihedral_nonhuman)=count_structures_groups(group_list['Nonhuman'])

            subList['All']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.dihedral==labelSelect,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
            subList['Human']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.specie=='Homo sapiens',Cluster.dihedral==labelSelect,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()
            subList['Nonhuman']=Cluster.query.filter(Cluster.group.contains(groupSelect),Cluster.specie!='Homo sapiens',Cluster.dihedral==labelSelect,Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).all()

            reprStr[tabs,'Any',labelSelect]=min_atom_missing(subList[tabs])
            nglList[tabs,reprStr[tabs,'Any',labelSelect]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'Any',labelSelect]).all()
            if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'Any',labelSelect]).first():
                dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'Any',labelSelect]).first().dfgnum)

            if labelSelect in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans'):
                pymolSession[tabs]=f'downloads/pymolSessions/{tabs}_{groupSelect}_DFGin_{labelSelect}_{ligand_type}.pse.zip'
                pymolScript[tabs]=f'downloads/pymolSessionScripts/{tabs}_{groupSelect}_DFGin_{labelSelect}_{ligand_type}.zip'
                pymolSessionRe[tabs]=f'downloads/pymolSessions/Repr_{tabs}_{groupSelect}_DFGin_{labelSelect}_{ligand_type}.pse.zip'
                pymolScriptRe[tabs]=f'downloads/pymolSessionScripts/Repr_{tabs}_{groupSelect}_DFGin_{labelSelect}_{ligand_type}.zip'
                csvFile[tabs]=f'downloads/text-files/{tabs}_{groupSelect}_DFGin_{labelSelect}_{ligand_type}.tab'
                write_text_file(group_list[tabs],csvFile[tabs])
                coordinateFiles[tabs]=f'downloads/coordinateFiles/{tabs}_{groupSelect}_DFGin_{labelSelect}_{ligand_type}'
            if labelSelect=='BABtrans':
                pymolSession[tabs]=f'downloads/pymolSessions/{tabs}_{groupSelect}_DFGinter_{labelSelect}_{ligand_type}.pse.zip'
                pymolScript[tabs]=f'downloads/pymolSessionScripts/{tabs}_{groupSelect}_DFGinter_{labelSelect}_{ligand_type}.zip'
                pymolSessionRe[tabs]=f'downloads/pymolSessions/Repr_{tabs}_{groupSelect}_DFGinter_{labelSelect}_{ligand_type}.pse.zip'
                pymolScriptRe[tabs]=f'downloads/pymolSessionScripts/Repr_{tabs}_{groupSelect}_DFGinter_{labelSelect}_{ligand_type}.zip'
                csvFile[tabs]=f'downloads/text-files/{tabs}_{groupSelect}_DFGinter_{labelSelect}_{ligand_type}.tab'
                write_text_file(group_list[tabs],csvFile[tabs])
                coordinateFiles[tabs]=f'downloads/coordinateFiles/{tabs}_{groupSelect}_DFGinter_{labelSelect}_{ligand_type}'
            if labelSelect=='BBAminus':
                pymolSession[tabs]=f'downloads/pymolSessions/{tabs}_{groupSelect}_DFGout_{labelSelect}_{ligand_type}.pse.zip'
                pymolScript[tabs]=f'downloads/pymolSessionScripts/{tabs}_{groupSelect}_DFGout_{labelSelect}_{ligand_type}.zip'
                pymolSessionRe[tabs]=f'downloads/pymolSessions/Repr_{tabs}_{groupSelect}_DFGout_{labelSelect}_{ligand_type}.pse.zip'
                pymolScriptRe[tabs]=f'downloads/pymolSessionScripts/Repr_{tabs}_{groupSelect}_DFGout_{labelSelect}_{ligand_type}.zip'
                csvFile[tabs]=f'downloads/text-files/{tabs}_{groupSelect}_DFGout_{labelSelect}_{ligand_type}.tab'
                write_text_file(group_list[tabs],csvFile[tabs])
                coordinateFiles[tabs]=f'downloads/coordinateFiles/{tabs}_{groupSelect}_DFGout_{labelSelect}_{ligand_type}'

        return render_template('Group-Dihedral.html',group_list=group_list,total_count=total_count,strCount_all=strCount_all,strCount_human=strCount_human,\
        strCount_nonhuman=strCount_nonhuman,reprStr=reprStr,geneCount_all=geneCount_all,geneCount_human=geneCount_human,geneCount_nonhuman=geneCount_nonhuman,\
        groupname=groupSelect,label=labelSelect,pymolSession=pymolSession,pymolScript=pymolScript,pymolSessionRe=pymolSessionRe,pymolScriptRe=pymolScriptRe,\
        clusterColor=clusterColor,coordinateFiles=coordinateFiles,csvFile=csvFile,nglList=nglList,dfgNumReprStr=dfgNumReprStr,ligand_type=ligand_type)


    if labelSelect in ('DFGinNone','DFGinterNone','DFGoutNone','NoneNone'):       #only for dihedral None

        for tabs in ('Human','All','Nonhuman'):
            group_list['All']=Cluster.query.filter(Cluster.spatial==labelSelect[0:-4],Cluster.dihedral=='None').order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['All']=Cluster.query.filter(Cluster.spatial==labelSelect[0:-4],Cluster.dihedral=='None').count()
            (strCount_all,geneCount_all,totalGroup_all,totalDihedral_all)=count_structures_groups(group_list['All'])

            group_list['Human']=Cluster.query.filter(Cluster.spatial==labelSelect[0:-4],Cluster.dihedral=='None',Cluster.specie=='Homo sapiens').order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Human']=Cluster.query.filter(Cluster.spatial==labelSelect[0:-4],Cluster.dihedral=='None',Cluster.specie=='Homo sapiens').count()
            (strCount_human,geneCount_human,totalGroup_human,totalDihedral_human)=count_structures_groups(group_list['Human'])

            group_list['Nonhuman']=Cluster.query.filter(Cluster.spatial==labelSelect[0:-4],Cluster.dihedral=='None',Cluster.specie!='Homo sapiens').order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Nonhuman']=Cluster.query.filter(Cluster.spatial==labelSelect[0:-4],Cluster.dihedral=='None',Cluster.specie!='Homo sapiens').count()
            (strCount_nonhuman,geneCount_nonhuman,totalGroup_nonhuman,totalDihedral_nonhuman)=count_structures_groups(group_list['Nonhuman'])

            subList['All']=Cluster.query.filter(Cluster.spatial==labelSelect[0:-4],Cluster.dihedral=='None').all()
            subList['Human']=Cluster.query.filter(Cluster.specie=='Homo sapiens',Cluster.spatial==labelSelect[0:-4],Cluster.dihedral=='None').all()
            subList['Nonhuman']=Cluster.query.filter(Cluster.specie!='Homo sapiens',Cluster.spatial==labelSelect[0:-4],Cluster.dihedral=='None').all()

            reprStr[tabs,labelSelect[0:-4],'None']=min_atom_missing(subList[tabs])
            nglList[tabs,reprStr[tabs,labelSelect[0:-4],'None']]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect[0:-4],'None']).all()
            if Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect[0:-4],'None']).first():
                dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,labelSelect[0:-4],'None']).first().dfgnum)

            pymolSession[tabs]=f'downloads/pymolSessions/{tabs}_Allgroups_'+labelSelect[0:-4]+f'_None_{ligand_type}.pse.zip'
            pymolScript[tabs]=f'downloads/pymolSessionScripts/{tabs}_Allgroups_'+labelSelect[0:-4]+f'_None_{ligand_type}.pse.zip'
            pymolSessionRe[tabs]=f'downloads/pymolSessions/Repr_{tabs}_Allgroups_'+labelSelect[0:-4]+f'_None_{ligand_type}.pse.zip'
            pymolScriptRe[tabs]=f'downloads/pymolSessionScripts/Repr_{tabs}_Allgroups_'+labelSelect[0:-4]+f'_None_{ligand_type}.pse.zip'

            csvFile[tabs]=f'downloads/text-files/{tabs}_Allgroups_'+labelSelect[0:-4]+f'_None_{ligand_type}.tab'
            write_text_file(group_list[tabs],csvFile[tabs])
            coordinateFiles[tabs]=f'downloads/coordinateFiles/{tabs}_Allgroups_'+labelSelect[0:-4]+f'_None_{ligand_type}'

        label=labelSelect[0:-4]

        return render_template('Spatial-NA.html',group_list=group_list,total_count=total_count,strCount_all=strCount_all,strCount_human=strCount_human,strCount_nonhuman=strCount_nonhuman,\
        reprStr=reprStr,nglList=nglList,dfgNumReprStr=dfgNumReprStr,geneCount_all=geneCount_all,geneCount_human=geneCount_human,geneCount_nonhuman=geneCount_nonhuman,label=label,pymolSession=pymolSession,\
        pymolScript=pymolScript,pymolSessionRe=pymolSessionRe,pymolScriptRe=pymolScriptRe,coordinateFiles=coordinateFiles,clusterColor=clusterColor,csvFile=csvFile)



@app.route('/<settings>/<queryname>')
def uniqueQuery(settings,queryname):
    gene_list=dict();total_count=dict();first_entry=dict();reprStr=dict();
    group_name=dict();protein_name=dict();uniprotid=dict();pymolSession=dict();pymolScript=dict();pymolSessionRe=dict();pymolScriptRe=dict();coordinateFiles=dict()
    nglList=dict();dfgNumReprStr=dict();subList=dict();dfgNumReprStr['Human']=list();dfgNumReprStr['All']=list();dfgNumReprStr['Nonhuman']=list()
    csvFile=dict();csvFile1=dict();csvFile2=dict()
    clusterColor=create_color_lists()

    if settings=='PDB':
        queryname=queryname.upper()
        if len(queryname)==5:
            queryname=queryname[0:-1]           #Remove chain so that only PDB id is always searched
        pdb_list=Cluster.query.filter(Cluster.pdb.contains(queryname)).all()

        formattedSeq='';nglList=dict()
        fhandleSeq=open(f'{pwd}/formattedSeq/'+queryname.upper()+'.html','r')
        formattedSeq=Markup(fhandleSeq.readline())
        for items in pdb_list:
            pdbReso=items.resolution;pdbGene=items.gene;pdbProtein=items.protein_name;pdbUniprot=items.uniprotid;pdbGroup=items.group
            pdbMutation=items.chain_mut
            pdbPhos=items.chain_phos
            pdbPseudo=items.status
            pdb_organism=items.specie
            domain_begin=items.domain_begin;domain_end=items.domain_end
            pheNum=int(items.dfgnum);aspNum=pheNum-1;xdfgNum=aspNum-1;glyNum=pheNum+1
            nglList[items.pdb]=Cluster.query.filter(Cluster.pdb==items.pdb).all()

        pymolSession=f'downloads/pymolSessions/{pdbGroup}_{pdbGene}_{queryname}.pse.zip'
        pymolScript=f'downloads/pymolSessionScripts/{pdbGroup}_{pdbGene}_{queryname}.zip'
        coordinateFiles=f'downloads/coordinateFiles/{pdbGroup}_{pdbGene}_{queryname}'

        return render_template('pdbs.html',queryname=queryname,pdb_list=pdb_list,pdbReso=pdbReso,pdbGene=pdbGene,pdbProtein=pdbProtein,\
        pdbUniprot=pdbUniprot,pdbGroup=pdbGroup,formattedSeq=formattedSeq,pdbMutation=pdbMutation,pdbPhos=pdbPhos,xdfgNum=xdfgNum,glyNum=glyNum,\
        clusterColor=clusterColor,pdbPseudo=pdbPseudo,domain_begin=domain_begin,domain_end=domain_end,pymolSession=pymolSession,pymolScript=pymolScript,\
        coordinateFiles=coordinateFiles,nglList=nglList,pheNum=pheNum,pdb_organism=pdb_organism)

    if settings=='GROUP':
        queryname=queryname.upper()
        return redirect(url_for('multipleQuery',groupSelect=queryname,labelSelect='All',ligTypeSelect='All'))

    if settings=='GENE' and queryname.upper() not in ('JAK1','JAK2','JAK3','EIF2AK4','RPS6KA1','RPS6KA3','RPS6KA5','RPS6KA6','TYK2'):
        queryname=queryname.upper()

        for tabs in ('Human','All','Nonhuman'):
                gene_list['Human']=Cluster.query.filter(func.lower(Cluster.gene)==func.lower(queryname),Cluster.specie=='Homo sapiens').order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
                total_count['Human']=Cluster.query.filter(func.lower(Cluster.gene)==func.lower(queryname),Cluster.specie=='Homo sapiens').count()
                first_entry['Human']=Cluster.query.filter(func.lower(Cluster.gene)==func.lower(queryname),Cluster.specie=='Homo sapiens').first()

                (strCount_human,geneCount_human,totalGroup_human,totalDihedral_human)=count_structures_groups(gene_list['Human'])

                gene_list['All']=Cluster.query.filter(func.lower(Cluster.gene)==func.lower(queryname)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
                total_count['All']=Cluster.query.filter(func.lower(Cluster.gene)==func.lower(queryname)).count()
                first_entry['All']=Cluster.query.filter(func.lower(Cluster.gene)==func.lower(queryname)).first()
                (strCount_all,geneCount_all,totalGroup_all,totalDihedral_all)=count_structures_groups(gene_list['All'])

                gene_list['Nonhuman']=Cluster.query.filter(func.lower(Cluster.gene)==func.lower(queryname),Cluster.specie!='Homo sapiens').order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
                total_count['Nonhuman']=Cluster.query.filter(func.lower(Cluster.gene)==func.lower(queryname),Cluster.specie!='Homo sapiens').count()
                first_entry['Nonhuman']=Cluster.query.filter(func.lower(Cluster.gene)==func.lower(queryname),Cluster.specie!='Homo sapiens').first()
                (strCount_nonhuman,geneCount_nonhuman,totalGroup_nonhuman,totalDihedral_nonhuman)=count_structures_groups(gene_list['Nonhuman'])

                if total_count['Human']>=1:      #try to assign groupname from humans, if not then other specie
                    geneName=first_entry['Human'].gene
                    group_name=first_entry['Human'].group;
                    protein_name=first_entry['Human'].protein_name;
                    uniprotid['Human']=first_entry['Human'].uniprotid
                    synonym=first_entry['Human'].synonym
                else:
                    geneName=first_entry['Nonhuman'].gene
                    group_name=first_entry['Nonhuman'].group;
                    protein_name=first_entry['Nonhuman'].protein_name;
                    uniprotid['Nonhuman']=first_entry['Nonhuman'].uniprotid

                for dihedral in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','None'):
                    subList['Human']=Cluster.query.filter(Cluster.gene==queryname,Cluster.spatial=='DFGin',Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens').all()
                    subList['All']=Cluster.query.filter(Cluster.gene==queryname,Cluster.spatial=='DFGin',Cluster.dihedral==dihedral).all()
                    subList['Nonhuman']=Cluster.query.filter(Cluster.gene==queryname,Cluster.spatial=='DFGin',Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens').all()

                    reprStr[tabs,'DFGin',dihedral]=min_atom_missing(subList[tabs])
                    nglList[tabs,reprStr[tabs,'DFGin',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGin',dihedral]).all()
                    if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGin',dihedral]).first():
                        dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGin',dihedral]).first().dfgnum)

                for dihedral in ('BABtrans','None'):
                    subList['Human']=Cluster.query.filter(Cluster.gene==queryname,Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens').all()
                    subList['All']=Cluster.query.filter(Cluster.gene==queryname,Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral).all()
                    subList['Nonhuman']=Cluster.query.filter(Cluster.gene==queryname,Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens').all()

                    reprStr[tabs,'DFGinter',dihedral]=min_atom_missing(subList[tabs])
                    nglList[tabs,reprStr[tabs,'DFGinter',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGinter',dihedral]).all()
                    if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGinter',dihedral]).first():
                        dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGinter',dihedral]).first().dfgnum)

                for dihedral in ('BBAminus','None'):
                    subList['Human']=Cluster.query.filter(Cluster.gene==queryname,Cluster.spatial=='DFGout',Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens').all()
                    subList['All']=Cluster.query.filter(Cluster.gene==queryname,Cluster.spatial=='DFGout',Cluster.dihedral==dihedral).all()
                    subList['Nonhuman']=Cluster.query.filter(Cluster.gene==queryname,Cluster.spatial=='DFGout',Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens').all()

                    reprStr[tabs,'DFGout',dihedral]=min_atom_missing(subList[tabs])
                    nglList[tabs,reprStr[tabs,'DFGout',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGout',dihedral]).all()
                    if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGout',dihedral]).first():
                        dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGout',dihedral]).first().dfgnum)

                subList['Human']=Cluster.query.filter(Cluster.gene==queryname,Cluster.spatial=='None',Cluster.specie=='Homo sapiens').all()
                subList['All']=Cluster.query.filter(Cluster.gene==queryname,Cluster.spatial=='None').all()
                subList['Nonhuman']=Cluster.query.filter(Cluster.gene==queryname,Cluster.spatial=='None',Cluster.specie!='Homo sapiens').all()

                reprStr[tabs,'None','None']=min_atom_missing(subList[tabs])
                nglList[tabs,reprStr[tabs,'None','None']]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).all()
                if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).first():
                    dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).first().dfgnum)

                pymolSession[tabs]=(f'downloads/pymolSessions/{tabs}_{group_name}_{queryname}.pse.zip')
                pymolScript[tabs]=(f'downloads/pymolSessionScripts/{tabs}_{group_name}_{queryname}.zip')
                pymolSessionRe[tabs]=(f'downloads/pymolSessions/Repr_{tabs}_{group_name}_{queryname}.pse.zip')
                pymolScriptRe[tabs]=(f'downloads/pymolSessionScripts/Repr_{tabs}_{group_name}_{queryname}.zip')

                csvFile[tabs]=f'downloads/text-files/{tabs}_{group_name}_{queryname}.tab'
                write_text_file(gene_list[tabs],csvFile[tabs])
                coordinateFiles[tabs]=f'downloads/coordinateFiles/{tabs}_{group_name}_{queryname}'

        return render_template('gene.html',first_entry=first_entry,gene_list=gene_list,geneName=geneName,protein_name=protein_name,group_name=group_name,\
        strCount_human=strCount_human,strCount_all=strCount_all,\
        strCount_nonhuman=strCount_nonhuman,reprStr=reprStr,pymolSession=pymolSession,pymolScript=pymolScript,pymolSessionRe=pymolSessionRe,pymolScriptRe=pymolScriptRe,coordinateFiles=coordinateFiles,clusterColor=clusterColor,\
        nglList=nglList,dfgNumReprStr=dfgNumReprStr,total_count=total_count,csvFile=csvFile)


    if settings=='UNIPROTID' and queryname.upper() not in ('JAK1_HUMAN','JAK2_HUMAN','JAK3_HUMAN','E2AK4_HUMAN','KS6A1_HUMAN','KS6A3_HUMAN',\
    'KS6A5_HUMAN','KS6A6_HUMAN','TYK2_HUMAN'):
        queryname=queryname.upper()
        first_entry=Cluster.query.filter(func.lower(Cluster.uniprotid)==func.lower(queryname)).first()
        geneName=first_entry.gene
        return redirect(url_for('uniqueQuery',settings='GENE',queryname=geneName))

    if settings=='GENE' and queryname.upper() in ('JAK1','JAK2','JAK3','EIF2AK4','RPS6KA1','RPS6KA3','RPS6KA5','RPS6KA6','TYK2'):
        queryname=queryname.upper()
        domain_1=(queryname+'_1');domain_2=(queryname+'_2')
        gene_list_1=dict();gene_list_2=dict();total_count_1=dict();total_count_2=dict();first_entry_1=dict();first_entry_2=dict();reprStr_1=dict();reprStr_2=dict();
        geneName='';group_name_1=dict();group_name_2=dict();uniprotid=dict();protein_name=dict();pymolSession_1=dict();pymolScript_1=dict();pymolSession_2=dict();\
        pymolScript_2=dict();pymolSessionRe_1=dict();pymolScriptRe_1=dict();pymolSessionRe_2=dict();\
        pymolScriptRe_2=dict();coordinateFiles=dict();nglList_1=dict();nglList_2=dict();dfgNumReprStr_1=dict();dfgNumReprStr_2=dict();subList=dict();\
        dfgNumReprStr_1['Human']=list();dfgNumReprStr_1['All']=list();dfgNumReprStr_1['Nonhuman']=list()
        dfgNumReprStr_2['Human']=list();dfgNumReprStr_2['All']=list();dfgNumReprStr_2['Nonhuman']=list();coordinateFiles_1=dict();coordinateFiles_2=dict()
        csvFile_1=dict();csvFile_2=dict()
        strCount_1_human='';geneCount_1_human='';totalGroup_1_human='';totalDihedral_1_human='';strCount_1_all='';geneCount_1_all='';totalGroup_1_all='';totalDihedral_1_all='';
        strCount_1_nonhuman='';geneCount_1_nonhuman='';totalGroup_1_nonhuman='';totalDihedral_1_nonhuman='';strCount_2_all='';geneCount_2_all='';totalGroup_2_all='';
        totalDihedral_2_all='';strCount_2_all='';geneCount_2_all='';totalGroup_2_all='';totalDihedral_2_all='';strCount_2_nonhuman='';geneCount_2_nonhuman='';
        totalGroup_2_nonhuman='';totalDihedral_2_nonhuman='';

        for tabs in ('Human','All','Nonhuman'):
            gene_list_1['Human']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_1),Cluster.specie=='Homo sapiens').order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            gene_list_2['Human']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_2),Cluster.specie=='Homo sapiens').order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count_1['Human']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_1),Cluster.specie=='Homo sapiens').count()
            total_count_2['Human']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_2),Cluster.specie=='Homo sapiens').count()

            gene_list_1['All']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_1)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            gene_list_2['All']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_2)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count_1['All']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_1)).count()
            total_count_2['All']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_2)).count()

            gene_list_1['Nonhuman']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_1),Cluster.specie!='Homo sapiens').order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            gene_list_2['Nonhuman']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_2),Cluster.specie!='Homo sapiens').order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count_1['Nonhuman']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_1),Cluster.specie!='Homo sapiens').count()
            total_count_2['Nonhuman']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_2),Cluster.specie!='Homo sapiens').count()

            if len(gene_list_1[tabs])>0:
                first_entry_1['Human']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_1),Cluster.specie=='Homo sapiens').first()
                first_entry_1['All']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_1)).first()
                first_entry_1['Nonhuman']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_1),Cluster.specie!='Homo sapiens').first()

                geneName=first_entry_1[tabs].gene;uniprotid[tabs]=first_entry_1[tabs].uniprotid       #These variable are the same in _1 and _2
                group_name_1[tabs]=first_entry_1[tabs].group;protein_name=first_entry_1[tabs].protein_name;synonym=first_entry_1[tabs].synonym
                (strCount_1_human,geneCount_1_human,totalGroup_1_human,totalDihedral_1_human)=count_structures_groups(gene_list_1['Human'])
                (strCount_1_all,geneCount_1_all,totalGroup_1_all,totalDihedral_1_all)=count_structures_groups(gene_list_1['All'])
                (strCount_1_nonhuman,geneCount_1_nonhuman,totalGroup_1_nonhuman,totalDihedral_1_nonhuman)=count_structures_groups(gene_list_1['Nonhuman'])

                for dihedral in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','None'):
                    subList['Human']=Cluster.query.filter(Cluster.domain==domain_1,Cluster.spatial=='DFGin',Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens').all()
                    subList['All']=Cluster.query.filter(Cluster.domain==domain_1,Cluster.spatial=='DFGin',Cluster.dihedral==dihedral).all()
                    subList['Nonhuman']=Cluster.query.filter(Cluster.domain==domain_1,Cluster.spatial=='DFGin',Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens').all()
                    reprStr_1[tabs,'DFGin',dihedral]=min_atom_missing(subList[tabs])
                    nglList_1[tabs,reprStr_1[tabs,'DFGin',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr_1[tabs,'DFGin',dihedral]).all()

                    if Cluster.query.filter(Cluster.pdb==reprStr_1[tabs,'DFGin',dihedral]).first():
                        dfgNumReprStr_1[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr_1[tabs,'DFGin',dihedral]).first().dfgnum)
                for dihedral in ('BABtrans','None'):
                    subList['Human']=Cluster.query.filter(Cluster.domain==domain_1,Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens').all()
                    subList['All']=Cluster.query.filter(Cluster.domain==domain_1,Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral).all()
                    subList['Nonhuman']=Cluster.query.filter(Cluster.domain==domain_1,Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens').all()
                    reprStr_1[tabs,'DFGinter',dihedral]=min_atom_missing(subList[tabs])
                    nglList_1[tabs,reprStr_1[tabs,'DFGinter',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr_1[tabs,'DFGinter',dihedral]).all()

                    if Cluster.query.filter(Cluster.pdb==reprStr_1[tabs,'DFGinter',dihedral]).first():
                        dfgNumReprStr_1[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr_1[tabs,'DFGinter',dihedral]).first().dfgnum)
                for dihedral in ('BBAminus','None'):
                    subList['Human']=Cluster.query.filter(Cluster.domain==domain_1,Cluster.spatial=='DFGout',Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens').all()
                    subList['All']=Cluster.query.filter(Cluster.domain==domain_1,Cluster.spatial=='DFGout',Cluster.dihedral==dihedral).all()
                    subList['Nonhuman']=Cluster.query.filter(Cluster.domain==domain_1,Cluster.spatial=='DFGout',Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens').all()
                    reprStr_1[tabs,'DFGout',dihedral]=min_atom_missing(subList[tabs])
                    nglList_1[tabs,reprStr_1[tabs,'DFGout',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr_1[tabs,'DFGout',dihedral]).all()

                    if Cluster.query.filter(Cluster.pdb==reprStr_1[tabs,'DFGout',dihedral]).first():
                        dfgNumReprStr_1[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr_1[tabs,'DFGout',dihedral]).first().dfgnum)

                subList['Human']=Cluster.query.filter(Cluster.domain==domain_1,Cluster.spatial=='None',Cluster.specie=='Homo sapiens').all()
                subList['All']=Cluster.query.filter(Cluster.domain==domain_1,Cluster.spatial=='None').all()
                subList['Nonhuman']=Cluster.query.filter(Cluster.domain==domain_1,Cluster.spatial=='None',Cluster.specie!='Homo sapiens').all()
                reprStr_1[tabs,'None','None']=min_atom_missing(subList[tabs])
                nglList_1[tabs,reprStr_1[tabs,'None','None']]=Cluster.query.filter(Cluster.pdb==reprStr_1[tabs,'None','None']).all()

                if Cluster.query.filter(Cluster.pdb==reprStr_1[tabs,'None','None']).first():
                    dfgNumReprStr_1[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr_1[tabs,'None','None']).first().dfgnum)

                pymolSession_1[tabs]=(f'downloads/pymolSessions/{tabs}_{group_name_1[tabs]}_{domain_1}.pse.zip')
                pymolScript_1[tabs]=(f'downloads/pymolSessionScripts/{tabs}_{group_name_1[tabs]}_{domain_1}.zip')
                pymolSessionRe_1[tabs]=(f'downloads/pymolSessions/Repr_{tabs}_{group_name_1[tabs]}_{domain_1}.pse.zip')
                pymolScriptRe_1[tabs]=(f'downloads/pymolSessionScripts/Repr_{tabs}_{group_name_1[tabs]}_{domain_1}.zip')
                csvFile_1[tabs]=f'downloads/text-files/{tabs}_{group_name_1[tabs]}_{domain_1}.tab'
                write_text_file(gene_list_1[tabs],csvFile_1[tabs])

                coordinateFiles_1[tabs]=f'downloads/coordinateFiles/{tabs}_{group_name_1[tabs]}_{domain_1}'

            if len(gene_list_2[tabs])>0:
                first_entry_2['Human']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_2),Cluster.specie=='Homo sapiens').first()
                first_entry_2['All']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_2)).first()
                first_entry_2['Nonhuman']=Cluster.query.filter(func.lower(Cluster.domain)==func.lower(domain_2),Cluster.specie!='Homo sapiens').first()

                geneName=first_entry_2[tabs].gene;uniprotid[tabs]=first_entry_2[tabs].uniprotid       #These variable are the same in _1 and _2
                group_name_2[tabs]=first_entry_2[tabs].group;protein_name=first_entry_2[tabs].protein_name
                (strCount_2_human,geneCount_2_human,totalGroup_2_human,totalDihedral_2_human)=count_structures_groups(gene_list_2['Human'])
                (strCount_2_all,geneCount_2_all,totalGroup_2_all,totalDihedral_2_all)=count_structures_groups(gene_list_2['All'])
                (strCount_2_nonhuman,geneCount_2_nonhuman,totalGroup_2_nonhuman,totalDihedral_2_nonhuman)=count_structures_groups(gene_list_2['Nonhuman'])

                for dihedral in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','None'):
                    subList['Human']=Cluster.query.filter(Cluster.domain==domain_2,Cluster.spatial=='DFGin',Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens').all()
                    subList['All']=Cluster.query.filter(Cluster.domain==domain_2,Cluster.spatial=='DFGin',Cluster.dihedral==dihedral).all()
                    subList['Nonhuman']=Cluster.query.filter(Cluster.domain==domain_2,Cluster.spatial=='DFGin',Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens').all()
                    reprStr_2[tabs,'DFGin',dihedral]=min_atom_missing(subList[tabs])
                    nglList_2[tabs,reprStr_2[tabs,'DFGin',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr_2[tabs,'DFGin',dihedral]).all()

                    if Cluster.query.filter(Cluster.pdb==reprStr_2[tabs,'DFGin',dihedral]).first():
                        dfgNumReprStr_2[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr_2[tabs,'DFGin',dihedral]).first().dfgnum)

                for dihedral in ('BABtrans','None'):
                    subList['Human']=Cluster.query.filter(Cluster.domain==domain_2,Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens').all()
                    subList['All']=Cluster.query.filter(Cluster.domain==domain_2,Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral).all()
                    subList['Nonhuman']=Cluster.query.filter(Cluster.domain==domain_2,Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens').all()
                    reprStr_2[tabs,'DFGinter',dihedral]=min_atom_missing(subList[tabs])
                    nglList_2[tabs,reprStr_2[tabs,'DFGinter',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr_2[tabs,'DFGinter',dihedral]).all()

                    if Cluster.query.filter(Cluster.pdb==reprStr_2[tabs,'DFGinter',dihedral]).first():
                        dfgNumReprStr_2[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr_2[tabs,'DFGinter',dihedral]).first().dfgnum)

                for dihedral in ('BBAminus','None'):
                    subList['Human']=Cluster.query.filter(Cluster.domain==domain_2,Cluster.spatial=='DFGout',Cluster.dihedral==dihedral,Cluster.specie=='Homo sapiens').all()
                    subList['All']=Cluster.query.filter(Cluster.domain==domain_2,Cluster.spatial=='DFGout',Cluster.dihedral==dihedral).all()
                    subList['Nonhuman']=Cluster.query.filter(Cluster.domain==domain_2,Cluster.spatial=='DFGout',Cluster.dihedral==dihedral,Cluster.specie!='Homo sapiens').all()
                    reprStr_2[tabs,'DFGout',dihedral]=min_atom_missing(subList[tabs])
                    nglList_2[tabs,reprStr_2[tabs,'DFGout',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr_2[tabs,'DFGout',dihedral]).all()

                    if Cluster.query.filter(Cluster.pdb==reprStr_2[tabs,'DFGout',dihedral]).first():
                        dfgNumReprStr_2[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr_2[tabs,'DFGout',dihedral]).first().dfgnum)

                subList['Human']=Cluster.query.filter(Cluster.domain==domain_2,Cluster.spatial=='None',Cluster.specie=='Homo sapiens').all()
                subList['All']=Cluster.query.filter(Cluster.domain==domain_2,Cluster.spatial=='None').all()
                subList['Nonhuman']=Cluster.query.filter(Cluster.domain==domain_2,Cluster.spatial=='None',Cluster.specie!='Homo sapiens').all()
                reprStr_2[tabs,'None','None']=min_atom_missing(subList[tabs])
                nglList_2[tabs,reprStr_2[tabs,'None','None']]=Cluster.query.filter(Cluster.pdb==reprStr_2[tabs,'None','None']).all()

                if Cluster.query.filter(Cluster.pdb==reprStr_2[tabs,'None','None']).first():
                    dfgNumReprStr_2[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr_2[tabs,'None','None']).first().dfgnum)

                pymolSession_2[tabs]=(f'downloads/pymolSessions/{tabs}_{group_name_2[tabs]}_{domain_2}.pse.zip')
                pymolScript_2[tabs]=(f'downloads/pymolSessionScripts/{tabs}_{group_name_2[tabs]}_{domain_2}.zip')
                pymolSessionRe_2[tabs]=(f'downloads/pymolSessions/Repr_{tabs}_{group_name_2[tabs]}_{domain_2}.pse.zip')
                pymolScriptRe_2[tabs]=(f'downloads/pymolSessionScripts/Repr_{tabs}_{group_name_2[tabs]}_{domain_2}.zip')
                csvFile_2[tabs]=f'downloads/text-files/{tabs}_{group_name_2[tabs]}_{domain_2}.tab'
                write_text_file(gene_list_2[tabs],csvFile_2[tabs])
                coordinateFiles_2[tabs]=f'downloads/coordinateFiles/{tabs}_{group_name_2[tabs]}_{domain_2}'

        return render_template('domain.html',first_entry_1=first_entry_1,first_entry_2=first_entry_2,geneName=geneName,protein_name=protein_name,\
        uniprotid=uniprotid,gene_list_1=gene_list_1,reprStr_1=reprStr_1,strCount_1_human=strCount_1_human,strCount_1_all=strCount_1_all,strCount_1_nonhuman=strCount_1_nonhuman,\
        total_count_1=total_count_1,pymolSession_1=pymolSession_1,pymolScript_1=pymolScript_1,gene_list_2=gene_list_2,reprStr_2=reprStr_2,strCount_2_human=strCount_2_human,\
        strCount_2_all=strCount_2_all,strCount_2_nonhuman=strCount_2_nonhuman,total_count_2=total_count_2,pymolSession_2=pymolSession_2,pymolScript_2=pymolScript_2,\
        coordinateFiles_1=coordinateFiles_1,coordinateFiles_2=coordinateFiles_2,clusterColor=clusterColor,csvFile_1=csvFile_1,csvFile_2=csvFile_2,nglList_1=nglList_1,\
        nglList_2=nglList_2,dfgNumReprStr_1=dfgNumReprStr_1,dfgNumReprStr_2=dfgNumReprStr_2)


    if settings=='UNIPROTID' and queryname.upper() in ('JAK1_HUMAN','JAK2_HUMAN','JAK3_HUMAN','E2AK4_HUMAN','KS6A1_HUMAN','KS6A3_HUMAN',\
    'KS6A5_HUMAN','KS6A6_HUMAN','TYK2_HUMAN'):
        queryname=queryname.upper()
        first_entry=Cluster.query.filter(func.lower(Cluster.uniprotid)==func.lower(queryname)).first()
        geneName=first_entry.gene
        return redirect(url_for('uniqueQuery',settings='GENE',queryname=geneName))


    if settings=='SPATIAL':
        return redirect(url_for('multipleQuery',groupSelect='All',labelSelect=queryname,ligTypeSelect='All'))

    if settings=='DIHEDRAL':
        return redirect(url_for('multipleQuery',groupSelect='All',labelSelect=queryname,ligTypeSelect='All'))

    if settings=='LIGAND':
        if queryname!='No_ligand':
            queryname=queryname.upper()
        ligand_list=dict();total_count=dict();geneList=dict();geneList['Human']=list();geneList['All']=list();geneList['Nonhuman']=list();
        reprStr=dict();nglList=dict();subList=dict();pymolSession=dict();pymolSession=dict();coordinateFiles=dict()
        dfgNumReprStr=dict();dfgNumReprStr['Human']=list();dfgNumReprStr['All']=list();dfgNumReprStr['Nonhuman']=list()

        for tabs in ('Human','All','Nonhuman'):
            ligand_list['Human']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.specie=='Homo sapiens').order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Human']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.specie=='Homo sapiens').count()
            (strCount_human,geneCount_human)=count_structures_all(ligand_list['Human'])

            ligand_list['All']=Cluster.query.filter(Cluster.ligand.contains(queryname)).order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['All']=Cluster.query.filter(Cluster.ligand.contains(queryname)).count()
            (strCount_all,geneCount_all)=count_structures_all(ligand_list['All'])

            ligand_list['Nonhuman']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.specie!='Homo sapiens').order_by(Cluster.specie).order_by(Cluster.group).order_by(Cluster.gene).order_by(Cluster.spatial).order_by(Cluster.dihedral).order_by(asc(Cluster.ligand_type)).all()
            total_count['Nonhuman']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.specie!='Homo sapiens').count()
            (strCount_nonhuman,geneCount_nonhuman)=count_structures_all(ligand_list['Nonhuman'])

            for spatial in ('DFGin','DFGinter','DFGout','None'):
                if spatial=='DFGin':
                    for dihedral in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','None'):
                        subList['Human']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.specie=='Homo sapiens',Cluster.spatial=='DFGin',Cluster.dihedral==dihedral).all()
                        subList['All']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.spatial=='DFGin',Cluster.dihedral==dihedral).all()
                        subList['Nonhuman']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.specie!='Homo sapiens',Cluster.spatial=='DFGin',Cluster.dihedral==dihedral).all()

                        reprStr[tabs,'DFGin',dihedral]=min_atom_missing(subList[tabs])
                        nglList[tabs,reprStr[tabs,'DFGin',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGin',dihedral]).all()
                        if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGin',dihedral]).first():
                            dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGin',dihedral]).first().dfgnum)
                if spatial=='DFGinter':
                    for dihedral in ('BABtrans','None'):
                        subList['Human']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.specie=='Homo sapiens',Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral).all()
                        subList['All']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral).all()
                        subList['Nonhuman']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.specie!='Homo sapiens',Cluster.spatial=='DFGinter',Cluster.dihedral==dihedral).all()

                        reprStr[tabs,'DFGinter',dihedral]=min_atom_missing(subList[tabs])
                        nglList[tabs,reprStr[tabs,'DFGinter',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGinter',dihedral]).all()
                        if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGinter',dihedral]).first():
                            dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGinter',dihedral]).first().dfgnum)
                if spatial=='DFGout':
                    for dihedral in ('BBAminus','None'):
                        subList['Human']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.specie=='Homo sapiens',Cluster.spatial=='DFGout',Cluster.dihedral==dihedral).all()
                        subList['All']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.spatial=='DFGout',Cluster.dihedral==dihedral).all()
                        subList['Nonhuman']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.specie!='Homo sapiens',Cluster.spatial=='DFGout',Cluster.dihedral==dihedral).all()

                        reprStr[tabs,'DFGout',dihedral]=min_atom_missing(subList[tabs])
                        nglList[tabs,reprStr[tabs,'DFGout',dihedral]]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGout',dihedral]).all()
                        if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGout',dihedral]).first():
                            dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'DFGout',dihedral]).first().dfgnum)

                subList['Human']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.specie=='Homo sapiens',Cluster.spatial=='None').all()
                subList['All']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.spatial=='None').all()
                subList['Nonhuman']=Cluster.query.filter(Cluster.ligand.contains(queryname),Cluster.specie!='Homo sapiens',Cluster.spatial=='None').all()

                reprStr[tabs,'None','None']=min_atom_missing(subList[tabs])
                nglList[tabs,reprStr[tabs,'None','None']]=Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).all()
                if Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).first():
                    dfgNumReprStr[tabs].append(Cluster.query.filter(Cluster.pdb==reprStr[tabs,'None','None']).first().dfgnum)

            for item in ligand_list[tabs]:
                geneList[tabs].append(item.gene)

            geneCount=len(set(geneList[tabs]))

            pymolSession[tabs]=(f'downloads/pymolSessions/{tabs}_{queryname}.pse.zip')
            pymolScript[tabs]=(f'downloads/pymolSessionScripts/{tabs}_{queryname}.zip')

            csvFile[tabs]=f'downloads/text-files/{tabs}_{queryname}.tab'
            write_text_file(ligand_list[tabs],csvFile[tabs])
            coordinateFiles[tabs]=f'downloads/coordinateFiles/{tabs}_{queryname}'

        return render_template('ligand.html',ligand_list=ligand_list,ligandname=queryname,geneCount_human=geneCount_human,geneCount_all=geneCount_all,\
        geneCount_nonhuman=geneCount_nonhuman,pymolSession=pymolSession,pymolScript=pymolScript,strCount_human=strCount_human,strCount_all=strCount_all,\
        strCount_nonhuman=strCount_nonhuman,reprStr=reprStr,total_count=total_count,clusterColor=clusterColor,coordinateFiles=coordinateFiles,\
        csvFile=csvFile,nglList=nglList,dfgNumReprStr=dfgNumReprStr)


@app.route('/webserver', methods=['GET','POST'])
def webserver():
    if request.method=='POST':
        if request.files:
            userpdb=request.files['userpdb.pdb']            #userpdb.pdb is the name in html file which is linked here, but if multiple users are accessing the server at the same time then the same filename can create a problem
            print(userpdb)
            if userpdb.filename=='':
                return 'File should have a valid name'
            if not allowed_filename(userpdb.filename):
                return 'File extension is not valid'
            # if '.cif' in userpdb.filename.lower():
            #     newFilename=userpdb.filename+'.cif'
            #     print(newFilename)
            # if '.pdb' in userpdb.filename.lower():
            #     newFilename=userpdb.filename+str(random.randint(1000,9999))+'.pdb'
            #     print(newFilename)
            #userpdb.save(os.path.join(UPLOAD_FOLDER,userpdb.filename))
            userpdb.save(os.path.join(UPLOAD_FOLDER,userpdb.filename))

            (group,chain_list,xdfg,dfg_asp,dfg_phe,xdfg_res,dfg_asp_res,dfg_phe_res,xdfg_phi,xdfg_psi,dfg_asp_phi,dfg_asp_psi,dfg_phe_phi,dfg_phe_psi,dfg_phe_chi1,chelix_conf,dfg_label,dfg_bkbone)=identify_state(pwd,userpdb.filename)
            print(userpdb.filename,group,dfg_label,dfg_bkbone)
            return render_template('conformation.html',userfilename=userpdb.filename,group=group,chain_list=chain_list,xdfg=xdfg,dfg_asp=dfg_asp,dfg_phe=dfg_phe,xdfg_res=xdfg_res,dfg_asp_res=dfg_asp_res,\
            dfg_phe_res=dfg_phe_res,xdfg_phi=xdfg_phi,xdfg_psi=xdfg_psi,dfg_asp_phi=dfg_asp_phi,dfg_asp_psi=dfg_asp_psi,\
            dfg_phe_phi=dfg_phe_phi,dfg_phe_psi=dfg_phe_psi,dfg_phe_chi1=dfg_phe_chi1,chelix_conf=chelix_conf,dfg_label=dfg_label,dfg_bkbone=dfg_bkbone)

    return render_template('webserver.html')


@app.route('/FDA')
def FDA():
    return render_template('FDA.html')

@app.route('/dunbrackLab')
def dunbrackLab():
    return redirect("http://dunbrack.fccc.edu")

@app.context_processor
def update_date():      #To update data in the footer; this route could be accessed by any html page
    fhandle_update_date=open(f'{pwd}/update-date.txt','r')
    read_date=fhandle_update_date.readline()
    return dict(updateDate=read_date)



###################################################
if __name__ == '__main__':
    app.run(debug=True)
