"""Generates code for EmJetHistoMaker
Conventions:

user_XXX functions are used to define specfic instances of elements.
generate_XXX functions are called from cog to generate code fragments.
All other functions are used internally."""

from histo_functions import *
import array

standalone = False
if __name__=='__main__':
    print "Running standalone mode - Printing to screen"
    standalone = True
else:
    try: import cog
    except: ImportError

def outputline(line):
    if not standalone:
        cog.outl(line)
    else:
        print(line)

############################################################
# Object definitions
############################################################
# namedtuple objects can be used as lightweight classes - similar to passive structs in C++

def user_define_bins():
    import math
    """Define bins to be used for variable width histograms in this function"""
    binning_dict = OrderedDict()
    binedges = map(lambda x: 10**x, compute_fixed_bins(100, -8., +2.))
    name = 'track_ipXY'          ; binning_dict[name] = VBins(len(binedges)-1, binedges, name)
    binedges = map(lambda x: 10**x, compute_fixed_bins(100, -2., +3.))
    name = 'track_ipSig'         ; binning_dict[name] = VBins(len(binedges)-1, binedges, name)
    binedges = map(lambda x: 10**x, compute_fixed_bins(200, -2., +8.))
    name = 'track_3dSig'         ; binning_dict[name] = VBins(len(binedges)-1, binedges, name)
    binedges = map(lambda x: 10**x, compute_fixed_bins(100, -3., +2.))
    name = 'track_distanceToJet' ; binning_dict[name] = VBins(len(binedges)-1, binedges, name)
    binedges = map(lambda x: 10**x, compute_fixed_bins(20, -3., +2.))
    name = 'vertex_Lxy'          ; binning_dict[name] = VBins(len(binedges)-1, binedges, name)
    name = 'vertex_mass'         ; binning_dict[name] = VBins(len(binedges)-1, binedges, name)
    # binedges = map(lambda x: 10**x, compute_fixed_bins(1000, -3+5., +2+5.))
    # binedges.insert(0, 0) # Add 0 as lower limit for first bin
    # name = 'sys_track_ipXY'      ; binning_dict[name] = VBins(len(binedges)-1, binedges, name)
    return binning_dict

def user_define_histos():
    """Define histograms in this function"""
    vbins = user_define_bins()
    histo_dict = OrderedDict()
    name = 'cutflow'                   ; histo_dict[name] = Histo1F(name , Bins( 25 , 0   ,  25  ) );
    name = 'cutflow2'                  ; histo_dict[name] = Histo1F(name , Bins( 25 , 0   ,  25  ) );
    name = 'weight_trigger'            ; histo_dict[name] = Histo1F(name , Bins(100 , 0   ,  10. ) );
    name = 'pdfshift'                  ; histo_dict[name] = Histo1F(name , Bins(100 , -1. ,  1.  ) );
    name = 'met_pt'                    ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 1000 ) )
    name = 'nJet'                      ; histo_dict[name] = Histo1F(name , Bins( 25 , 0   ,  25  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'nEmerging'                 ; histo_dict[name] = Histo1F(name , Bins( 25 , 0   ,  25  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'nJet'                      ; histo_dict[name] = Histo1F(name , Bins( 25 , 0   ,  25  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'pileup_nVtx'               ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 100  ) );
    name = 'nVtx'                      ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 100  ) );
    name = 'ht'                        ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 2500 ) )
    name = 'ht4'                       ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 2500 ) )
    name = 'norm_ht4'                  ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 2500 ) )
    name = 'PFHT900'                   ; histo_dict[name] = Histo1F(name , Bins(  2 , 0   ,    2 ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'PFHT900__PFHT475'          ; histo_dict[name] = Histo1F(name , Bins(  2 , 0   ,    2 ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'sigmaPt'                   ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 1500 ) )
    name = 'sigmaPt2'                  ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 1500 ) )
    name = 'deltaPt'                   ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 1000 ) )
    name = 'vertex_Lxy'                ; histo_dict[name] = Histo1F(name , vbins[name]             )
    name = 'vertex_mass'               ; histo_dict[name] = Histo1F(name , vbins[name]             )
    name = 'jet_N'                     ; histo_dict[name] = Histo1F(name , Bins( 25 , 0   ,  25  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'jet_pt'                    ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 1000 ) )
    name = 'jet_eta'                   ; histo_dict[name] = Histo1F(name , Bins(100 , -5  , 5    ) )
    name = 'jet_phi'                   ; histo_dict[name] = Histo1F(name , Bins(100 , -5  , 5    ) )
    name = 'jet_nTrack'                ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) )
    name = 'jet_medianIP'              ; histo_dict[name] = Histo1F(name , vbins['track_ipXY']     )
    name = 'jet_maxIP'                 ; histo_dict[name] = Histo1F(name , vbins['track_ipXY']     )
    name = 'jet_nTrackPostCut'         ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) )
    name = 'jet_medianIPPostCut'       ; histo_dict[name] = Histo1F(name , vbins['track_ipXY']     )
    name = 'jet_medianLogIpSig'        ; histo_dict[name] = Histo1F(name , Bins(100 , -5  , 5    ) )
    name = 'jet_medLogIpSig'           ; histo_dict[name] = Histo1F(name , Bins(100 , -5  , 5    ) )
    name = 'jet_alphaMax'              ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz'           ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz1um'        ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz2um'        ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz5um'        ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz10um'       ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz20um'       ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz50um'       ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz100um'      ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz200um'      ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz500um'      ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz1mm'        ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz2mm'        ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz5mm'        ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz1cm'        ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz2cm'        ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_alphaMax_dz5cm'        ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_nDarkPions'            ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'jet_nVertex'               ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'jet_prompt_frac_E'         ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 2.   ) )
    name = 'jet_prompt_frac'           ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 2.   ) )
    name = 'jet_disp_frac'             ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 2.   ) )
    name = 'jet_alphaMaxNeg'           ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_cef'                   ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_nef'                   ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_chf'                   ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_nhf'                   ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_phf'                   ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_missInnerHit_frac'     ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'sumMedianLogIpSig'         ; histo_dict[name] = Histo1F(name , Bins(100 , -25 ,  25  ) )
    # name = 'sys_track_ipXY'            ; histo_dict[name] = Histo1D(name , vbins['sys_track_ipXY'] )
    name = 'sys_log_track_ipXY'        ; histo_dict[name] = Histo1D(name , Bins(2000, -8. , 2.   ) )
    name = 'sys_log_track_ipXYSig'     ; histo_dict[name] = Histo1D(name , Bins(2000, -5. , 5.   ) )
    name = 'sys_track_3dSig'           ; histo_dict[name] = Histo1D(name , Bins(2000, -1. , 50.  ) )
    name = 'systest_log_track_ipXY'    ; histo_dict[name] = Histo1F(name , Bins(2000, -8. , 2.   ) )
    name = 'systest_track_3dSig'       ; histo_dict[name] = Histo1D(name , Bins(2000, -1. , 50.  ) )
    name = 'test_jet_medianIP'         ; histo_dict[name] = Histo1F(name , vbins['track_ipXY']     )
    name = 'test_jet_medianAbsIP'      ; histo_dict[name] = Histo1F(name , vbins['track_ipXY']     )
    name = 'track_pt'                  ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 10.  ) )
    name = 'track_eta'                 ; histo_dict[name] = Histo1F(name , Bins(100 , -5  , 5    ) )
    name = 'track_phi'                 ; histo_dict[name] = Histo1F(name , Bins(100 , -5  , 5    ) )
    name = 'track_ipXY'                ; histo_dict[name] = Histo1F(name , vbins['track_ipXY']     )
    name = 'track_ipXYb'               ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 2.  ) )
    name = 'track_ipSig'               ; histo_dict[name] = Histo1F(name , vbins['track_ipSig']    )
    name = 'track_logIpSig'            ; histo_dict[name] = Histo1F(name , Bins(100 , -5  ,   5  ) )
    name = 'track_ip3DSig'             ; histo_dict[name] = Histo1F(name , vbins['track_ipSig']    )
    name = 'track_tk3DSig'             ; histo_dict[name] = Histo1F(name , vbins['track_ipSig']    )
    name = 'track_nHits'               ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_nMissInnerHits'      ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_missHitFrac'         ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  ,   5. ) )
    name = 'track_nTrkLayers'          ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_nMissInnerTrkLayers' ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_nPxlLayers'          ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_nMissInnerPxlLayers' ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_nNetMissInnerLayers' ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_missLayerFrac'       ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  ,   5. ) )
    name = 'track_quality'             ; histo_dict[name] = Histo1F(name , Bins(200 , 0.  , 200  ) )
    name = 'track_algo'                ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) )
    name = 'track_originalAlgo'        ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) )
    name = 'track_dRToJetAxis'         ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 0.5  ) )
    name = 'track_distanceToJet'       ; histo_dict[name] = Histo1F(name , vbins[name]             )
    name = 'vertex_Lxy'                ; histo_dict[name] = Histo1F(name , vbins['vertex_Lxy']     )
    name = 'vertex_mass'               ; histo_dict[name] = Histo1F(name , vbins['vertex_mass']    )
    name = 'vertex_chi2'               ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) )
    name = 'vertex_ndof'               ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  ,  20  ) )
    name = 'vertex_pt2sum'             ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) )
    name = 'vertex_nTracks'            ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) )
    name = 'overlap_jet_minDeltaR'     ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 10.  ) )
    name = 'overlap_jet_nDupTracks'    ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) )
    name = 'metUnc_met_delta'          ; histo_dict[name] = Histo1F(name , Bins(100 ,-500 , 500  ) )
    name = 'metUnc_jet_met_pt'         ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 1000 ) )
    name = 'metUnc_jet_met_delta_up'   ; histo_dict[name] = Histo1F(name , Bins(100 , -20 , 20   ) )
    name = 'metUnc_jet_met_delta_dn'   ; histo_dict[name] = Histo1F(name , Bins(100 , -20 , 20   ) )
    # name = 'metUnc_delta_met_pt'       ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 1000 ) )

    # Define 2D histograms from ordered pairs of 1D histograms
    histo_2d_dict = OrderedDict()
    h = histo_combine1Dto2D( histo_dict['track_eta']   , histo_dict['track_pt']    , ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['track_eta']   , histo_dict['track_phi']   , ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['track_phi']  , histo_dict['track_ipXY']   , ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['track_nHits'] , histo_dict['track_ipXY'] , ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['track_nHits'] , histo_dict['track_ipSig'] , ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['jet_alphaMax'], histo_dict['jet_medLogIpSig'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['jet_alphaMax'], histo_dict['jet_prompt_frac'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['jet_alphaMax'], histo_dict['jet_disp_frac'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['jet_alphaMax'], histo_dict['jet_maxIP'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['jet_pt'], histo_dict['jet_alphaMax'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['jet_pt'], histo_dict['jet_disp_frac'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['nVtx'], histo_dict['jet_alphaMax'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['ht'], histo_dict['jet_alphaMax'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['ht'], histo_dict['PFHT900'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['ht'], histo_dict['PFHT900__PFHT475'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['ht4'], histo_dict['PFHT900'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['ht4'], histo_dict['PFHT900__PFHT475'], ); histo_2d_dict[h.name] = h
    histo_dict.update(histo_2d_dict)

    # Track plot variations
    histo_clone_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if name[:6]=='track_' or name=='nTrack':
            histo_clone = clone_object(histo, postfix='TKprompt')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='TKdisplaced')
            histo_clone_dict[histo_clone.name] = histo_clone
    histo_dict.update(histo_clone_dict)

    # Jet plot variations
    histo_clone_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if name[:4]=='jet_' or name[:6]=='track_' or name[:7]=='vertex_' or name=='nJet':
            histo_clone = clone_object(histo, postfix='JTbasic')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='JTalphaMax')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='JTipcut')
            histo_clone_dict[histo_clone.name] = histo_clone
    histo_dict.update(histo_clone_dict)

    # Define signal histos, by copying all jet, track and vertex histograms
    histo_sig_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if name[:4]=='jet_' or name[:6]=='track_' or name[:7]=='vertex_':
            histo_sig = clone_object(histo, postfix='sig')
            histo_sig_dict[histo_sig.name] = histo_sig
    histo_dict.update(histo_sig_dict)

    # Event plot variations
    histo_clone_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if name[:4]=='jet_' or name[:6]=='track_' or name[:7]=='vertex_':
            histo_clone = clone_object(histo, postfix='EVTpvpass')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='EVTkinematic')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='EVTallpass')
            histo_clone_dict[histo_clone.name] = histo_clone
    histo_dict.update(histo_clone_dict)

    # Systematic variations
    histo_clone_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if name.startswith('cutflow') or name.startswith('test'):
            histo_clone = clone_object(histo, postfix='PdfUp')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='PdfDn')
            histo_clone_dict[histo_clone.name] = histo_clone
    histo_dict.update(histo_clone_dict)
    histo_clone_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if name.startswith('cutflow') or name.startswith('test'):
            histo_clone = clone_object(histo, postfix='TriggerUp')
            histo_clone_dict[histo_clone.name] = histo_clone
    histo_dict.update(histo_clone_dict)
    histo_clone_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if name.startswith('cutflow') or name.startswith('test'):
            histo_clone = clone_object(histo, postfix='ModelingUp')
            histo_clone_dict[histo_clone.name] = histo_clone
    histo_dict.update(histo_clone_dict)
    histo_clone_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if name.startswith('cutflow') or name.startswith('test'):
            histo_clone = clone_object(histo, postfix='PileupUp')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='PileupDn')
            histo_clone_dict[histo_clone.name] = histo_clone
    histo_dict.update(histo_clone_dict)
    histo_clone_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if name.startswith('cutflow') or name.startswith('test'):
            histo_clone = clone_object(histo, postfix='PtUp')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='PtDn')
            histo_clone_dict[histo_clone.name] = histo_clone
    histo_dict.update(histo_clone_dict)

    return histo_dict


def user_define_histo_vectors():
    """Define histogram vectors in this function"""
    histo_dict = user_define_histos()
    histo_vector_dict = OrderedDict()
    return histo_vector_dict

def construct_bin_str(bins):
    # Fixed width binning
    if type(bins).__name__ == 'Bins':
        return ', {nBins}, {min}, {max}'.format(nBins=bins.nBins, min=bins.min, max=bins.max)
    # Variable width binning
    if type(bins).__name__ == 'VBins':
        return ', nBins_{binname}, bins_{binname}'.format(nBins=bins.nBins, binname=bins.binname)

def construct_histo_init(obj):
    name    = obj.name
    typestr = get_type_str(obj)
    binstr  = construct_bin_str(obj.binsX)
    if hasattr(obj, 'binsY'): binstr += construct_bin_str(obj.binsY)
    if hasattr(obj, 'binsZ'): binstr += construct_bin_str(obj.binsZ)
    output = '{name} = new {typestr}("{name}", "{name}" {binstr});'.format(name=name, typestr=typestr, binstr=binstr)
    return output

def construct_histo_map_init(obj):
    name    = obj.name
    typestr = get_type_str(obj)
    mapstr  = get_map_str(obj)
    binstr  = construct_bin_str(obj.binsX)
    if hasattr(obj, 'binsY'): binstr += construct_bin_str(obj.binsY)
    if hasattr(obj, 'binsZ'): binstr += construct_bin_str(obj.binsZ)
    output = '{mapstr}["{name}"] = new {typestr}("{name}", "{name}" {binstr});'.format(name=name, typestr=typestr, mapstr=mapstr, binstr=binstr)
    return output

def generate_vbin_decl():
    for name, bins in user_define_bins().iteritems():
        nBins=bins.nBins
        bin_array_str = str(bins.binedges)[1:-1] # Remove first and last character (square brackets)
        bin_array_str = '{ ' + bin_array_str + ' }'
        outputline( 'const int nBins_{binname} = {nBins}; const double bins_{binname}[{nBins}+1] = {binedges};'.format(nBins=bins.nBins, binname=bins.binname, binedges=bin_array_str) )

def generate_histo_map_init():
    for name, histo in user_define_histos().iteritems():
        outputline( construct_histo_map_init(histo) )

def calculate_index(histo_dict):
    """Takes histo_dict and calculates index in vectors based on their order in dictionary and their type (1D vs 2D)"""
    index_1d = 0
    index_2d = 0
    index_1d_double = 0
    histo_id_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if type(histo).__name__ == 'Histo1F':
            histo_id_dict[name] = index_1d
            index_1d += 1
        if type(histo).__name__ == 'Histo1D':
            histo_id_dict[name] = index_1d_double
            index_1d_double += 1
        elif type(histo).__name__ == 'Histo2F':
            histo_id_dict[name] = index_2d
            index_2d += 1
    return histo_id_dict

def generate_histo_index():
    histo_id_dict = calculate_index( user_define_histos() )
    for name, histo in user_define_histos().iteritems():
        template_str = 'if (name=="{name}") return {index};'
        index = histo_id_dict[name]
        output_str = template_str.format(name=name, index=index)
        outputline(output_str)

if __name__=='__main__':
    print 'generate_histo_index():'
    generate_histo_index()
    print''
    print 'generate_vbin_decl():'
    generate_vbin_decl()
    print''
    print 'generate_histo_map_init():'
    generate_histo_map_init()
    print''

