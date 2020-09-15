
from KinasesApp import db, Cluster
#from flask_sqlalchemy import SQLAlchemy
import sys
import pandas as pd
from datetime import datetime

filename=sys.argv[1]
#fhandle=open(filename,'r')
today=str(datetime.now())[0:10].strip()

#df=pd.read_csv(f'Kinases_df-{today}.csv',sep='\t')
df=pd.read_csv(filename,sep='\t')
db.drop_all()
db.create_all()
for i in df.index:
    if df.at[i,'Ligand']=='No_ligand':
        df.at[i,'Ligand_label']='No_ligand'

    pdb1=Cluster(pdb=df.at[i,'PDBid'],\
    protein_name=df.at[i,'Protein'],\
    uniprotid=df.at[i,'UniprotID'],\
    chain_length=df.at[i,'ChainLen'],\
    str_begin=df.at[i,'StrBegin'],\
    str_end=df.at[i,'StrEnd'],\
    specie=df.at[i,'Specie'],\
    resolution=df.at[i,'Resolution'],\
    method=df.at[i,'Method'],\
    rvalue=df.at[i,'Rvalue'],\
    freervalue=df.at[i,'FreeRvalue'],\
    uniseq=df.at[i,'UniSeq'],\
    gene=df.at[i,'Gene'],\
    domain=df.at[i,'Domain'],\
    group=df.at[i,'Group'],\
    uniprotacc=df.at[i,'UniAcc'],\
    domain_begin=df.at[i,'DomainBegin'],\
    domain_end=df.at[i,'DomainEnd'],\
    alkres=df.at[i,'ALKres'],\
    alknum=df.at[i,'ALKnum'],\
    rreres=df.at[i,'RREres'],\
    rrenum=df.at[i,'RREnum'],\
    hrdres=df.at[i,'HRDres'],\
    hrdnum=df.at[i,'HRDnum'],\
    dfgres=df.at[i,'DFGres'],\
    dfgnum=df.at[i,'DFGnum'],\
    aperes=df.at[i,'APEres'],\
    apenum=df.at[i,'APEnum'],\
    dfgaspres=df.at[i,'DFG_Aspres'],\
    gtknum=df.at[i,'GTKnum'],\
    gtkres=df.at[i,'GTKres'],\
    hinge1=df.at[i,'Hinge1'],\
    status=df.at[i,'Status'],\
    #x_n_edia=df.at[i,'X_N_Edia'],\
    x_o_edia=df.at[i,'X_O_Edia'],\
    #asp_n_edia=df.at[i,'Asp_N_Edia'],\
    asp_o_edia=df.at[i,'Asp_O_Edia'],\
    #phe_n_edia=df.at[i,'Phe_N_Edia'],\
    phe_o_edia=df.at[i,'Phe_O_Edia'],\
    #gly_n_edia=df.at[i,'Gly_N_Edia'],\
    gly_o_edia=df.at[i,'Gly_O_Edia'],\
    x_phi=df.at[i,'X_Phi'],\
    x_psi=df.at[i,'X_Psi'],\
    asp_phi=df.at[i,'Asp_Phi'],\
    asp_psi=df.at[i,'Asp_Psi'],\
    asp_chi1=df.at[i,'Asp_Chi1'],\
    asp_chi2=df.at[i,'Asp_Chi2'],\
    phe_phi=df.at[i,'Phe_Phi'],\
    phe_psi=df.at[i,'Phe_Psi'],\
    phe_chi1=df.at[i,'Phe_Chi1'],\
    phe_chi2=df.at[i,'Phe_Chi2'],\
    gly_phi=df.at[i,'Gly_Phi'],\
    gly_psi=df.at[i,'Gly_Psi'],\
    domainBreak=df.at[i,'DomainBreak'],\
    loopBreak=int(df.at[i,'LoopBreak']),\
    synonym='X',\
    x_mut=df.at[i,'X_mut'],\
    asp_mut=df.at[i,'Asp_mut'],\
    phe_mut=df.at[i,'Phe_mut'],\
    gly_mut=df.at[i,'Gly_mut'],\
    chain_mut=df.at[i,'Chain_mut'],\
    chain_phos=df.at[i,'Chain_phos'],\
    modified_aa=df.at[i,'Modified_aa'],\
    first_obs_res=df.at[i,'First_obs_res'],\
    last_obs_res=df.at[i,'Last_obs_res'],\
    xdfg_resolved=df.at[i,'XDFGresolved'],\
    ligand=df.at[i,'Ligand'],\
    strseq=df.at[i,'StrSeq'],\
    seqres=df.at[i,'SeqRes'],\
    #lig_rre4=df.at[i,'Lig_RRE4'],\
    #lig_hinge=df.at[i,'Lig_Hinge'],\
    lys_glu=df.at[i,'Lys_Glu'],\
    phe_glu4=df.at[i,'Phe_Glu4'],\
    phe_lys=df.at[i,'Phe_Lys'],\
    chelix=df.at[i,'C-helix'],\
    spatial=df.at[i,'Spatial'],\
    dihedral=df.at[i,'Dihedral'],\
    #dihedral_nochi1=df.at[i,'Dihedral_NoChi1'],\
    color=df.at[i,'Color'],\
    ligand_type=df.at[i,'Ligand_label'])


    db.session.add(pdb1)
    db.session.commit()
