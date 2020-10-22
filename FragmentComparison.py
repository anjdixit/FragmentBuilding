def GenerateConfs(infile = "", outfile = ""):

    
# Initialize OEOmega object and set the max configurations per molecule. Set strict stero and atom types to false to keep restrictions loose
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(20)
    omega.SetStrictStereo(False)
    omega.SetStrictAtomTypes(False)


    ofs = oechem.oemolostream('%s.oeb.gz' % outfile)


    istream = oechem.oemolistream('%s.sdf' % infile) 

    for oemol in istream.GetOEMols():
        omega(oemol)
        oechem.OEWriteMolecule(ofs, oemol) 
    istream.close()
    ofs.close()

def SingleOverlayComparison(Ref_smiles = "", infile = "", outfile=""): # takes in SMILES of reference molecule, input file and returns tanimoto list of comparisons.
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(20)
    omega.SetStrictStereo(False)
    omega.SetStrictAtomTypes(False)

#set output file
    ofs = oechem.oemolostream(outfile)
    ofs.SetFormat(oechem.OEFormat_SDF)
#create mol for reference molecule
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, Ref_smiles)
    omega(mol)
    oechem.OEWriteMolecule(ofs, mol) 
    tanimoto = []


# Setup ROCS to provide specified number of conformers per hit
# Add our reference molecule as the one we are fitting    

    istream = oechem.oemolistream(infile)

    for item in istream.GetOEMols():
        options = OEROCSOptions()
        options.SetNumBestHits(3)
        rocs = OEROCS(options)
        rocs.AddMolecule(item)
        # Loop over results and output
        for res in rocs.Overlay(mol):                #add in overlay "test" molecule
            outmol = res.GetOverlayConf()              #Use GetOverlayConf to get just the best; GetOverlayConfs for all
            oeshape.OERemoveColorAtoms(outmol)
            oechem.OEAddExplicitHydrogens(outmol)

            tanimoto.append(res.GetTanimotoCombo())
            oechem.OEAddSDData(outmol,"Tanimoto", str(res.GetTanimotoCombo()))
            oechem.OEWriteMolecule(ofs, outmol)       
    istream.close()
    ofs.close()
    return tanimoto
	
def MatrixOverlayComparison(infile = "", outfile = "") # takes in file with conformers to be used to generate an nxn matrix of overlays. Returns an nxn numpy array of reverse tanimoto scores (or distances).
	
	istream = oechem.oemolistream(infile)
	ofs = oechem.oemolostream(outfile)
	ofs.SetFormat(oechem.OEFormat_SDF)
	
	mols = []
	
	for mol in istream.GetOEMols():
		mols.append(oechem.OEMol(mol)
		
	overlay_dist = np.zeros(len(mols))
	
	for i, mol1 in enumerate(mols):
		options = OEROCSOptions()
		options.SetNumBestHits(1)
		#options.SetConfsPerHit(10)
		rocs = OEROCS(options)
		rocs.AddMolecule(oechem.OEMol(mol1))                             #add in reference molecule
		
		
		for j, mol2 in enumerate(mols):
			if i==j:                   #take out values on the diagonal, or all the "self" pairs
				continue
			# Loop over results and output
	
			for res in rocs.Overlay(mol2):                #add in overlay "test" molecule
				outmol = res.GetOverlayConf()              #Use GetOverlayConf to get just the best; GetOverlayConfs for all
				oeshape.OERemoveColorAtoms(outmol)
				oechem.OEAddExplicitHydrogens(outmol)
				outmol.SetTitle(oechem.OEMolToSmiles(outmol))
	   
				oechem.OEAddSDData(outmol,"Tanimoto", str(res.GetTanimotoCombo()))
				oechem.OEAddSDData(outmol, "Molecule 1",  oechem.OEMolToSmiles(mol1))
				oechem.OEAddSDData(outmol, "Molecule 2",  oechem.OEMolToSmiles(mol2))
				
				# overlay_Confs.append(oechem.OEMol(outmol))
				overlay_dist[i,j] = (2 - res.GetTanimotoCombo()) if res.GetTanimotoCombo() <= 2 else 0 #get "distance" from tanimoto score by reversing it; sill set values of the diagonals to 0 (same molecules are compared)
				
				print(overlay_dist[i,j])
				oechem.OEWriteMolecule(ofs, outmol)
				
            
	ofs.close()
	istream.close()
	return overlay_dist
