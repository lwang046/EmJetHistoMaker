"""Generates code for EmJetHistoMaker
Conventions:

user_XXX functions are used to define specfic instances of elements.
generate_XXX functions are called from cog to generate code fragments.
All other functions are used internally."""

from histo_functions import *

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
    binedges = map(lambda x: 10**x, compute_fixed_bins(100, -3., +2.))
    name = 'track_distanceToJet' ; binning_dict[name] = VBins(len(binedges)-1, binedges, name)
    binedges = map(lambda x: 10**x, compute_fixed_bins(20, -3., +2.))
    name = 'vertex_Lxy'          ; binning_dict[name] = VBins(len(binedges)-1, binedges, name)
    name = 'vertex_mass'         ; binning_dict[name] = VBins(len(binedges)-1, binedges, name)
    return binning_dict

def user_define_histos():
    """Define histograms in this function"""
    vbins = user_define_bins()
    histo_dict = OrderedDict()
    name = 'cutflow'                   ; histo_dict[name] = Histo1F(name , Bins( 25 , 0   ,  25  ) );
    name = 'nJet'                      ; histo_dict[name] = Histo1F(name , Bins( 25 , 0   ,  25  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'nJetPostCut'               ; histo_dict[name] = Histo1F(name , Bins( 25 , 0   ,  25  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'ht4'                       ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 2500 ) )
    name = 'ht'                        ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 2500 ) )
    name = 'sigmaPt'                   ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 1500 ) )
    name = 'sigmaPt2'                  ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 1500 ) )
    name = 'deltaPt'                   ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 1000 ) )
    name = 'vertex_Lxy'                ; histo_dict[name] = Histo1F(name , vbins[name]             )
    name = 'vertex_mass'               ; histo_dict[name] = Histo1F(name , vbins[name]             )
    name = 'jet_pt'                    ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 1000 ) )
    name = 'jet_eta'                   ; histo_dict[name] = Histo1F(name , Bins(100 , -5  , 5    ) )
    name = 'jet_nTracks'               ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'jet_medianLogIpSig'        ; histo_dict[name] = Histo1F(name , Bins(100 , -5  , 5    ) )
    name = 'jet_medLogIpSig'           ; histo_dict[name] = Histo1F(name , Bins(100 , -5  , 5    ) )
    name = 'jet_alphaMax'              ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_nDarkPions'            ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'jet_nVertex'               ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'jet_prompt_frac_E'         ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 2.   ) )
    name = 'jet_displaced_frac'        ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 2.   ) )
    name = 'jet_alphaMaxNeg'           ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_cef'                   ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'jet_nef'                   ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 1.   ) )
    name = 'sumMedianLogIpSig'         ; histo_dict[name] = Histo1F(name , Bins(100 , -25 ,  25  ) )
    name = 'track_pt'                  ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 10.  ) )
    name = 'track_eta'                 ; histo_dict[name] = Histo1F(name , Bins(100 , -5  , 5    ) )
    name = 'track_phi'                 ; histo_dict[name] = Histo1F(name , Bins(100 , -5  , 5    ) )
    name = 'track_ipXY'                ; histo_dict[name] = Histo1F(name , vbins['track_ipXY']     )
    name = 'track_ipSig'               ; histo_dict[name] = Histo1F(name , vbins['track_ipSig']    )
    name = 'track_logIpSig'            ; histo_dict[name] = Histo1F(name , Bins(100 , -5  ,   5  ) )
    name = 'track_nHits'               ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_nMissHits'           ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_missHitFrac'         ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  ,   5. ) )
    name = 'track_nTrkLayers'          ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_nMissInnerTrkLayers' ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_nPxlLayers'          ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_nMissInnerPxlLayers' ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_nNetMissInnerLayers' ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  , 100  ) ); histo_dict[name] = offset(histo_dict[name])
    name = 'track_missLayerFrac'       ; histo_dict[name] = Histo1F(name , Bins(100 , 0.  ,   5. ) )
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

    # Define 2D histograms from ordered pairs of 1D histograms
    histo_2d_dict = OrderedDict()
    h = histo_combine1Dto2D( histo_dict['track_eta']   , histo_dict['track_pt']    , ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['track_eta']   , histo_dict['track_phi']   , ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['track_phi']  , histo_dict['track_ipXY']   , ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['track_nHits'] , histo_dict['track_ipXY'] , ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['track_nHits'] , histo_dict['track_ipSig'] , ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['jet_alphaMax'], histo_dict['jet_medLogIpSig'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['jet_alphaMax'], histo_dict['jet_prompt_frac_E'], ); histo_2d_dict[h.name] = h
    h = histo_combine1Dto2D( histo_dict['jet_alphaMaxNeg'], histo_dict['jet_displaced_frac'], ); histo_2d_dict[h.name] = h
    histo_dict.update(histo_2d_dict)

    # Define signal histos, by copying all jet, track and vertex histograms
    histo_sig_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if name[:4]=='jet_' or name[:6]=='track_' or name[:7]=='vertex_':
            histo_sig = clone_object(histo, postfix='sig')
            histo_sig_dict[histo_sig.name] = histo_sig
    histo_dict.update(histo_sig_dict)

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

if __name__=='__main__':
    print 'generate_vbin_decl():'
    generate_vbin_decl()
    print''
    print 'generate_histo_map_init():'
    generate_histo_map_init()
    print''

