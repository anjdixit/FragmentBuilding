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

def SingleComparisonOverlay(Ref_smiles = "", infile = "", outfile=""): # takes in SMILES of reference molecule, input file basename (oeb.gz) and 
                                                                       # returns tanimoto list of comparisons.
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(20)
    omega.SetStrictStereo(False)
    omega.SetStrictAtomTypes(False)

#set output file
    ofs = oechem.oemolostream('%s.sdf' % outfile)
    ofs.SetFormat(oechem.OEFormat_SDF)
#create mol for reference molecule
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, Ref_smiles)
    omega(mol)
    oechem.OEWriteMolecule(ofs, mol) 
    tanimoto = []


# Setup ROCS to provide specified number of conformers per hit
# Add our reference molecule as the one we are fitting    

    istream = oechem.oemolistream('%s.oeb.gz'% infile)

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

