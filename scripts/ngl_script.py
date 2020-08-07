#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 16:07:55 2020

@author: vivekmodi
"""

def write_ngl_script(pwd,pdb_domain_dict,pdb_dfgnum_dict):
    for pdbs in pdb_domain_dict:
        domain_name=pdb_domain_dict[pdbs]
        filename=('/static/js/pdbs.js')
        pdbfile=(pwd+'kinasechains_renumber_uniprot/'+pdbs+'.cif.gz')
        fhandle_js=open(filename,'w')
        
        fhandle_js.write('var stage = new NGL.Stage("viewport");\n')
        fhandle_js.write('')
stage.loadFile("rcsb://1CRN").then(function (component) {
  // add a "cartoon" representation to the structure component
  component.addRepresentation("cartoon");
  // provide a "good" view of the structure
  component.autoView();
});
