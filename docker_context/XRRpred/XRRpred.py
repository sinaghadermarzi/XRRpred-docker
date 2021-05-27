#!/usr/bin/env python 
from residuelevel_scores import *
from prot_feature_implement import *
import os
from joblib import load
from save_output import *
import json
from gen_taskid import *






avg_dict = {'frac_maxc_sidechainclass_aliphatic': 0.372129110881029,
 'frac_maxc_sidechainclass_basic': 0.143291101097883,
 'frac_maxc_sidechainclass_amide': 0.082766338887211,
 'frac_maxc_sidechainclass_acid': 0.13114828510877802,
 'frac_maxc_sidechainclass_sulfurcontaining': 0.037406487867963,
 'frac_maxc_sidechainclass_aromatic': 0.115939311028702,
 'frac_maxc_sidechainpolarity_nonpolar': 0.412181109941026,
 'frac_maxc_sidechainpolarity_basicpolar': 0.143291101097883,
 'frac_maxc_sidechainpolarity_polar': 0.234451431319987,
 'frac_maxc_sidechainpolarity_acidicpolar': 0.13114828510877802,
 'frac_maxc_sidechaincharge_neutral': 0.7716965787775341,
 'frac_maxc_sidechaincharge_positive': 0.11271386444942402,
 'frac_maxc_sidechaincharge_negative': 0.13114828510877802,
 'frac_maxc_heavyHFRYW': 0.168500459873723,
 'frac_maxc_lightGASPV': 0.330258580629664,
 'frac_maxc_hydrophob': 0.653062848365934,
 'frac_maxc_tinyAGCS': 0.229474861950856,
 'frac_maxc_smallAGCSVTDNP': 0.4997351002772361,
 'frac_minc_sidechainclass_aliphatic': 0.3628806280937821,
 'frac_minc_sidechainclass_basic': 0.138315223354676,
 'frac_minc_sidechainclass_amide': 0.079875739414485,
 'frac_minc_sidechainclass_acid': 0.127498969936997,
 'frac_minc_sidechainclass_sulfurcontaining': 0.035405932408455,
 'frac_minc_sidechainclass_aromatic': 0.111495961184629,
 'frac_minc_sidechainpolarity_nonpolar': 0.399497418330541,
 'frac_minc_sidechainpolarity_basicpolar': 0.138315223354676,
 'frac_minc_sidechainpolarity_polar': 0.22701672681366894,
 'frac_minc_sidechainpolarity_acidicpolar': 0.127498969936997,
 'frac_minc_sidechaincharge_neutral': 0.7480987622334409,
 'frac_minc_sidechaincharge_positive': 0.109050500356687,
 'frac_minc_sidechaincharge_negative': 0.127498969936997,
 'frac_minc_heavyHFRYW': 0.162329754041055,
 'frac_minc_lightGASPV': 0.321265234695562,
 'frac_minc_hydrophob': 0.632878202882243,
 'frac_minc_tinyAGCS': 0.222881110119983,
 'frac_minc_smallAGCSVTDNP': 0.4851327043823879,
 'frac_avgc_sidechainclass_aliphatic': 0.367472931765584,
 'frac_avgc_sidechainclass_basic': 0.140776361920822,
 'frac_avgc_sidechainclass_amide': 0.081327573363982,
 'frac_avgc_sidechainclass_acid': 0.12932519589755698,
 'frac_avgc_sidechainclass_sulfurcontaining': 0.036391925243247004,
 'frac_avgc_sidechainclass_aromatic': 0.113698311569453,
 'frac_avgc_sidechainpolarity_nonpolar': 0.405766762623174,
 'frac_avgc_sidechainpolarity_basicpolar': 0.140776361920822,
 'frac_avgc_sidechainpolarity_polar': 0.230727663141293,
 'frac_avgc_sidechainpolarity_acidicpolar': 0.12932519589755698,
 'frac_avgc_sidechaincharge_neutral': 0.759815829837782,
 'frac_avgc_sidechaincharge_positive': 0.110858974264658,
 'frac_avgc_sidechaincharge_negative': 0.12932519589755698,
 'frac_avgc_heavyHFRYW': 0.165378558702358,
 'frac_avgc_lightGASPV': 0.3257144272357461,
 'frac_avgc_hydrophob': 0.6428801862692529,
 'frac_avgc_tinyAGCS': 0.226125485446072,
 'frac_avgc_smallAGCSVTDNP': 0.492369693970657,
 'no_unique_chains': 1.06216880633707,
 'sum_chain_length': 603.392242556679,
 'average_chain_length': 272.715674863317,
 'shortest_chain_length': 267.960885004097,
 'longest_chain_length': 277.664791040699,
 'num_chains': 2.33652007648184,
 'stdev_chains': 4.67646860355472,
 'frac_maxc_top3Flexbilityidx': 0.183177747647654,
 'frac_maxc_top3bvalueidx': 0.19130100071318,
 'frac_maxc_top3topidp': 0.177643013562146,
 'frac_maxc_top3foldunfold': 0.085362074380242,
 'frac_maxc_top3disprot': 0.15722767757852898,
 'frac_maxc_bottom3Flexbilityidx': 0.07577956575643201,
 'frac_maxc_bottom3bvalueidx': 0.065182845739159,
 'frac_maxc_bottom3topidp': 0.085362074380242,
 'frac_maxc_bottom3foldunfold': 0.17597094904891394,
 'frac_maxc_bottom3disprot': 0.083217376138485,
 'frac_minc_top3Flexbilityidx': 0.178073896473445,
 'frac_minc_top3bvalueidx': 0.185715677437256,
 'frac_minc_top3topidp': 0.17206577699134198,
 'frac_minc_top3foldunfold': 0.08223123818663901,
 'frac_minc_top3disprot': 0.15208857446106305,
 'frac_minc_bottom3Flexbilityidx': 0.072790949478617,
 'frac_minc_bottom3bvalueidx': 0.062271764566335,
 'frac_minc_bottom3topidp': 0.08223123818663901,
 'frac_minc_bottom3foldunfold': 0.170942464718314,
 'frac_minc_bottom3disprot': 0.080055924713942,
 'frac_avgc_top3Flexbilityidx': 0.180598367321488,
 'frac_avgc_top3bvalueidx': 0.188503923029308,
 'frac_avgc_top3topidp': 0.174839301297201,
 'frac_avgc_top3foldunfold': 0.08378092391328901,
 'frac_avgc_top3disprot': 0.15465246467554195,
 'frac_avgc_bottom3Flexbilityidx': 0.07426735285401201,
 'frac_avgc_bottom3bvalueidx': 0.063703073272822,
 'frac_avgc_bottom3topidp': 0.08378092391328901,
 'frac_avgc_bottom3foldunfold': 0.17343479177828194,
 'frac_avgc_bottom3disprot': 0.08161728843478701,
 'frac_maxc_AA_A': 0.080395075812336,
 'frac_maxc_AA_R': 0.052561148845019987,
 'frac_maxc_AA_N': 0.043028959266423986,
 'frac_maxc_AA_D': 0.059012337397318,
 'frac_maxc_AA_C': 0.013404883925345,
 'frac_maxc_AA_Q': 0.039737379620787,
 'frac_maxc_AA_E': 0.072135947711459,
 'frac_maxc_AA_G': 0.07160426140531599,
 'frac_maxc_AA_H': 0.03057723664846,
 'frac_maxc_AA_I': 0.056808720437448,
 'frac_maxc_AA_L': 0.09448680086805802,
 'frac_maxc_AA_K': 0.060152715604405,
 'frac_maxc_AA_M': 0.024001603942618,
 'frac_maxc_AA_F': 0.038774190038121,
 'frac_maxc_AA_P': 0.045354350246281,
 'frac_maxc_AA_S': 0.064070640807859,
 'frac_maxc_AA_T': 0.054030339058487,
 'frac_maxc_AA_W': 0.013003771775693002,
 'frac_maxc_AA_Y': 0.033584112566429,
 'frac_maxc_AA_V': 0.06883425235787201,
 'frac_minc_AA_A': 0.078247354218319,
 'frac_minc_AA_R': 0.050833792856426,
 'frac_minc_AA_N': 0.041636234444504,
 'frac_minc_AA_D': 0.05735226154721201,
 'frac_minc_AA_C': 0.012443373748087,
 'frac_minc_AA_Q': 0.038239504969981,
 'frac_minc_AA_E': 0.070146708389784,
 'frac_minc_AA_G': 0.06988784206980701,
 'frac_minc_AA_H': 0.029264722997990007,
 'frac_minc_AA_I': 0.055300400492049,
 'frac_minc_AA_L': 0.092319894091237,
 'frac_minc_AA_K': 0.058216707500262,
 'frac_minc_AA_M': 0.022962558660369,
 'frac_minc_AA_F': 0.037516240344441,
 'frac_minc_AA_P': 0.043702361101296,
 'frac_minc_AA_S': 0.062302540083769,
 'frac_minc_AA_T': 0.052435599947024014,
 'frac_minc_AA_W': 0.012312150473808,
 'frac_minc_AA_Y': 0.032402847368391004,
 'frac_minc_AA_V': 0.067125137222369,
 'frac_avgc_AA_A': 0.079290887170636,
 'frac_avgc_AA_R': 0.051680247132905,
 'frac_avgc_AA_N': 0.042335682853888,
 'frac_avgc_AA_D': 0.05818129332181,
 'frac_avgc_AA_C': 0.012913822831029,
 'frac_avgc_AA_Q': 0.038991890510093,
 'frac_avgc_AA_E': 0.071143902575747,
 'frac_avgc_AA_G': 0.070736826866774,
 'frac_avgc_AA_H': 0.029917387656164,
 'frac_avgc_AA_I': 0.05605510828002299,
 'frac_avgc_AA_L': 0.093404016417148,
 'frac_avgc_AA_K': 0.059178727131754,
 'frac_avgc_AA_M': 0.023478102412219,
 'frac_avgc_AA_F': 0.038140893118057,
 'frac_avgc_AA_P': 0.04451667158969999,
 'frac_avgc_AA_S': 0.063183948577632,
 'frac_avgc_AA_T': 0.053224467728182005,
 'frac_avgc_AA_W': 0.012648357323737,
 'frac_avgc_AA_Y': 0.032991673471496,
 'frac_avgc_AA_V': 0.067986093031004,
 'idx_maxc_maxswnum20_top3_flex': 7.4174269325320985,
 'idx_maxc_maxswnum20_top3_bvalue': 7.73280524446872,
 'idx_maxc_maxswnum20_top3_topidp': 7.410871346626611,
 'idx_maxc_maxswnum20_top3_foldunfold': 4.548320131111719,
 'idx_maxc_maxswnum20_top3_disprot': 6.809833378858238,
 'iul_maxc_avg': 0.246068701327007,
 'iul_maxc_maxswsum10': 0.43233964887884496,
 'iul_maxc_maxswsum20': 0.40684710708218896,
 'iul_maxc_maxswsum25': 0.389929835252498,
 'iul_maxc_maxswsum30': 0.373899706906844,
 'iul_maxc_maxswsum35': 0.360624695874988,
 'iul_maxc_frac1sbin': 0.06736660306685201,
 'iul_maxc_maxcons1sbin': 6.353947009013928,
 'iul_maxc_maxcons1sbinpl': 0.03971540160268,
 'iul_minc_avg': 0.23764178958663695,
 'iul_minc_maxswsum10': 0.422678503141217,
 'iul_minc_maxswsum20': 0.39773500961566,
 'iul_minc_maxswsum25': 0.3814604254742899,
 'iul_minc_maxswsum30': 0.365886041104312,
 'iul_minc_maxswsum35': 0.353072863733515,
 'iul_minc_frac1sbin': 0.057800527009396,
 'iul_minc_maxcons1sbin': 5.82403714832013,
 'iul_minc_maxcons1sbinpl': 0.032627468499886,
 'iul_avgc_avg': 0.241788017035438,
 'iul_avgc_maxswsum10': 0.42747601598832,
 'iul_avgc_maxswsum20': 0.402265799200194,
 'iul_avgc_maxswsum25': 0.385671085022774,
 'iul_avgc_maxswsum30': 0.369871377250949,
 'iul_avgc_maxswsum35': 0.35682982004499897,
 'iul_avgc_frac1sbin': 0.062396742040285,
 'iul_avgc_maxcons1sbin': 6.07894139344222,
 'iul_avgc_maxcons1sbinpl': 0.036034846573264,
 'ius_maxc_avg': 0.207141308772294,
 'ius_maxc_maxswsum10': 0.4082860122668929,
 'ius_maxc_maxswsum20': 0.388356700562849,
 'ius_maxc_maxswsum25': 0.372828193702099,
 'ius_maxc_maxswsum30': 0.358822977078159,
 'ius_maxc_maxswsum35': 0.34563556738178103,
 'ius_maxc_frac1sbin': 0.045858718086897,
 'ius_maxc_maxcons1sbin': 5.32515706091232,
 'ius_maxc_maxcons1sbinpl': 0.032262652084242,
 'ius_minc_avg': 0.19976118363477904,
 'ius_minc_maxswsum10': 0.399821658265253,
 'ius_minc_maxswsum20': 0.380220781830185,
 'ius_minc_maxswsum25': 0.365243524149284,
 'ius_minc_maxswsum30': 0.35163626159232503,
 'ius_minc_maxswsum35': 0.338850792418895,
 'ius_minc_frac1sbin': 0.037730323191156,
 'ius_minc_maxcons1sbin': 4.86189565692434,
 'ius_minc_maxcons1sbinpl': 0.026207452953715,
 'ius_avgc_avg': 0.203382080623271,
 'ius_avgc_maxswsum10': 0.404017944414838,
 'ius_avgc_maxswsum20': 0.384270199836538,
 'ius_avgc_maxswsum25': 0.369013326901622,
 'ius_avgc_maxswsum30': 0.355210478608204,
 'ius_avgc_maxswsum35': 0.3422213518456311,
 'ius_avgc_frac1sbin': 0.041626547283118986,
 'ius_avgc_maxcons1sbin': 5.08276593963565,
 'ius_avgc_maxcons1sbinpl': 0.02910745417353,
 'iul_maxc_avg_asa_th1': 0.25511519061793,
 'iul_maxc_maxswsum10_asa_th1': 0.42581674465023,
 'iul_maxc_maxswsum20_asa_th1': 0.405074529768064,
 'iul_maxc_maxswsum25_asa_th1': 0.39140088449743504,
 'iul_maxc_maxswsum30_asa_th1': 0.377535152517586,
 'iul_maxc_maxswsum35_asa_th1': 0.365460831383071,
 'iul_maxc_frac1sbin_asa_th1': 0.075174843322048,
 'iul_minc_avg_asa_th1': 0.246681505473145,
 'iul_minc_maxswsum10_asa_th1': 0.4161351006849199,
 'iul_minc_maxswsum20_asa_th1': 0.395942797808818,
 'iul_minc_maxswsum25_asa_th1': 0.382903817281708,
 'iul_minc_maxswsum30_asa_th1': 0.369477016370962,
 'iul_minc_maxswsum35_asa_th1': 0.357861123805928,
 'iul_minc_frac1sbin_asa_th1': 0.065232533667509,
 'iul_avgc_avg_asa_th1': 0.250835068546708,
 'iul_avgc_maxswsum10_asa_th1': 0.4209426071536479,
 'iul_avgc_maxswsum20_asa_th1': 0.400484891666846,
 'iul_avgc_maxswsum25_asa_th1': 0.387130635099224,
 'iul_avgc_maxswsum30_asa_th1': 0.37348435700904703,
 'iul_avgc_maxswsum35_asa_th1': 0.361638671333339,
 'iul_avgc_frac1sbin_asa_th1': 0.070014809288809,
 'ius_maxc_avg_asa_th1': 0.217515009062726,
 'ius_maxc_maxswsum10_asa_th1': 0.402769729428362,
 'ius_maxc_maxswsum20_asa_th1': 0.386643465416079,
 'ius_maxc_maxswsum25_asa_th1': 0.3736641575160311,
 'ius_maxc_maxswsum30_asa_th1': 0.361285484880149,
 'ius_maxc_maxswsum35_asa_th1': 0.349413319213422,
 'ius_maxc_frac1sbin_asa_th1': 0.052630634203375,
 'ius_minc_avg_asa_th1': 0.210091488879298,
 'ius_minc_maxswsum10_asa_th1': 0.3942534297988271,
 'ius_minc_maxswsum20_asa_th1': 0.378514299968606,
 'ius_minc_maxswsum25_asa_th1': 0.3660806785289461,
 'ius_minc_maxswsum30_asa_th1': 0.354070486259511,
 'ius_minc_maxswsum35_asa_th1': 0.34258990417832896,
 'ius_minc_frac1sbin_asa_th1': 0.044118556139339,
 'ius_avgc_avg_asa_th1': 0.213734119043914,
 'ius_avgc_maxswsum10_asa_th1': 0.3984725895668729,
 'ius_avgc_maxswsum20_asa_th1': 0.38256044372356496,
 'ius_avgc_maxswsum25_asa_th1': 0.369852671356427,
 'ius_avgc_maxswsum30_asa_th1': 0.35766022849711604,
 'ius_avgc_maxswsum35_asa_th1': 0.345982345358105,
 'ius_avgc_frac1sbin_asa_th1': 0.04820514187892,
 'iul_maxc_avg_asa_th2': 0.251219763744821,
 'iul_maxc_maxswsum10_asa_th2': 0.429799745494118,
 'iul_maxc_maxswsum20_asa_th2': 0.406356803911243,
 'iul_maxc_maxswsum25_asa_th2': 0.39114535264289096,
 'iul_maxc_maxswsum30_asa_th2': 0.376223018874064,
 'iul_maxc_maxswsum35_asa_th2': 0.363512077448495,
 'iul_maxc_frac1sbin_asa_th2': 0.071692269177246,
 'iul_minc_avg_asa_th2': 0.24279507939203104,
 'iul_minc_maxswsum10_asa_th2': 0.42013448044262,
 'iul_minc_maxswsum20_asa_th2': 0.397244750749153,
 'iul_minc_maxswsum25_asa_th2': 0.38266284868868,
 'iul_minc_maxswsum30_asa_th2': 0.368180122028664,
 'iul_minc_maxswsum35_asa_th2': 0.35593722993895704,
 'iul_minc_frac1sbin_asa_th2': 0.061939601102249,
 'iul_avgc_avg_asa_th2': 0.246941016270552,
 'iul_avgc_maxswsum10_asa_th2': 0.424935195675216,
 'iul_avgc_maxswsum20_asa_th2': 0.40177759048072703,
 'iul_avgc_maxswsum25_asa_th2': 0.38688266362024093,
 'iul_avgc_maxswsum30_asa_th2': 0.37218099114384295,
 'iul_avgc_maxswsum35_asa_th2': 0.359705299120404,
 'iul_avgc_frac1sbin_asa_th2': 0.06662803546528599,
 'ius_maxc_avg_asa_th2': 0.213094964795364,
 'ius_maxc_maxswsum10_asa_th2': 0.406220475682018,
 'ius_maxc_maxswsum20_asa_th2': 0.3879419918024461,
 'ius_maxc_maxswsum25_asa_th2': 0.373676398160084,
 'ius_maxc_maxswsum30_asa_th2': 0.360487341481002,
 'ius_maxc_maxswsum35_asa_th2': 0.347966728060615,
 'ius_maxc_frac1sbin_asa_th2': 0.04953779032144,
 'ius_minc_avg_asa_th2': 0.205704799097125,
 'ius_minc_maxswsum10_asa_th2': 0.3977472455069289,
 'ius_minc_maxswsum20_asa_th2': 0.379824350789489,
 'ius_minc_maxswsum25_asa_th2': 0.3660974806874779,
 'ius_minc_maxswsum30_asa_th2': 0.353288367508649,
 'ius_minc_maxswsum35_asa_th2': 0.34116478112771503,
 'ius_minc_frac1sbin_asa_th2': 0.041236000420202985,
 'ius_avgc_avg_asa_th2': 0.209329447361504,
 'ius_avgc_maxswsum10_asa_th2': 0.4019480056592221,
 'ius_avgc_maxswsum20_asa_th2': 0.3838662025606121,
 'ius_avgc_maxswsum25_asa_th2': 0.369866252967875,
 'ius_avgc_maxswsum30_asa_th2': 0.356869755159508,
 'ius_avgc_maxswsum35_asa_th2': 0.344544876637768,
 'ius_avgc_frac1sbin_asa_th2': 0.045218971686542,
 'iul_maxc_maxminswsum10': 0.075168682178232,
 'iul_maxc_maxminswsum20': 0.097446670267647,
 'iul_maxc_maxminswsum25': 0.108240358695053,
 'iul_maxc_maxminswsum30': 0.119785694559506,
 'iul_maxc_maxminswsum35': 0.129347386218749,
 'iul_minc_minswsum10': 0.388747937721935,
 'iul_minc_minswsum20': 0.320755904433118,
 'iul_minc_minswsum25': 0.289062331517851,
 'iul_minc_minswsum30': 0.2587596428242661,
 'iul_minc_minswsum35': 0.23420216515011,
 'ius_maxc_maxminswsum10': 0.042341039457669,
 'ius_maxc_maxminswsum20': 0.05585550675966599,
 'ius_maxc_maxminswsum25': 0.063244886818084,
 'ius_maxc_maxminswsum30': 0.07056847600465899,
 'ius_maxc_maxminswsum35': 0.077808005494268,
 'ius_minc_minswsum10': 0.399458235454795,
 'ius_minc_minswsum20': 0.346426461723928,
 'ius_minc_minswsum25': 0.319527468779227,
 'ius_minc_minswsum30': 0.295517297247538,
 'ius_minc_minswsum35': 0.27323358120896696,
 'iul_maxc_maxminswsum10_asa_th1': 0.065906734433732,
 'iul_maxc_maxminswsum20_asa_th1': 0.090914678845523,
 'iul_maxc_maxminswsum25_asa_th1': 0.103107128194091,
 'iul_maxc_maxminswsum30_asa_th1': 0.115650509929201,
 'iul_maxc_maxminswsum35_asa_th1': 0.126106504618605,
 'iul_minc_minswsum10_asa_th1': 0.391023669046968,
 'iul_minc_minswsum20_asa_th1': 0.32697365804549605,
 'iul_minc_minswsum25_asa_th1': 0.297582823672772,
 'iul_minc_minswsum30_asa_th1': 0.268285553977123,
 'iul_minc_minswsum35_asa_th1': 0.243799311708125,
 'ius_maxc_maxminswsum10_asa_th1': 0.036732664255224,
 'ius_maxc_maxminswsum20_asa_th1': 0.052547916857137,
 'ius_maxc_maxminswsum25_asa_th1': 0.060576643313396,
 'ius_maxc_maxminswsum30_asa_th1': 0.068421753069855,
 'ius_maxc_maxminswsum35_asa_th1': 0.076117201049807,
 'ius_minc_minswsum10_asa_th1': 0.400626092272474,
 'ius_minc_minswsum20_asa_th1': 0.350839183843404,
 'ius_minc_minswsum25_asa_th1': 0.32578876131564005,
 'ius_minc_minswsum30_asa_th1': 0.30252886555181296,
 'ius_minc_minswsum35_asa_th1': 0.28078305145500393}

























def calculate_features(prot):
    row  = {}
    for f in features_to_compute:
        row[f.__name__] = f(prot)
    return row

def normalize_by_scaler(in_df,features,scaler):
    X = in_df.loc[:,features].values
    X = scaler.transform(X)
    out_df = pandas.DataFrame(X,index = in_df.index, columns=features)
    return out_df

def predict(model,features_df,feature_set,min_val, max_val):
    test_X = features_df.loc[:,feature_set].values
    test_pred = model.predict(test_X)
    for j in range(len(test_pred)):
        if test_pred[j]> max_val:
            test_pred[j] = max_val
        elif test_pred[j]<min_val:
            test_pred[j] = min_val
    return list(test_pred)


this_script_directory = os.path.dirname(os.path.realpath(__file__))


in_path = sys.argv[1]
in_path = os.path.abspath(in_path)


# taskid = gen_taskid("run")
# with open(this_script_directory+"/"+taskid+".txt","w") as cmdfile:
#     cmdfile.writelines(in_path)


with open(in_path) as infile:
    input_text = infile.read()



# with open(this_script_directory+"/"+taskid+"_txt.txt","w") as cmdfile:
#     cmdfile.writelines(input_text)

chains = dict()
pids = []
entries = input_text[1:].rstrip("\n").split(">")
for entry in entries:
    lines = entry.split("\n")
    entry_id = lines[0]
    seq = "".join(lines[1:])
    if "_" in entry_id:
        spl = entry_id.split("_")
        pid = spl[0]
    else:
        pid = entry_id
    if pid in chains:
        chains[pid].append(seq)
    else:
        chains[pid]= [seq]
        pids.append(pid)

# taskid = gen_taskid("xrrpred")
# with open(this_script_directory+"/"+taskid+"_chains.txt","w") as cmdfile:
#     cmdfile.writelines(str(chains))

# with open(this_script_directory+"/"+taskid+"_pids.txt","w") as cmdfile:
#     cmdfile.writelines(str(pids))


profile_exit_code = []
profile_stop_chain= []
chain_data = []
for pid in pids:
    try:
        current_prot  = []
        for seq in chains[pid]:
            current_chain = {}
            current_chain["id"] = pid
            current_chain["sequence"] = seq
            current_chain["iupred_short"] = iupred_short(seq)
            current_chain["iupred_long"] = iupred_long(seq)
            current_chain["asaquick"] = asaquick(seq)
            current_prot.append(current_chain)
        profile_exit_code.append("success")
    except:
        if ("iupred_long" in current_chain):
            profile_exit_code.append("ASAquick")
        elif ("iupred_short" in current_chain):
            profile_exit_code.append("IUPred long version")
        else:
            profile_exit_code.append("IUPred short version")
    finally:
        profile_stop_chain.append(len(current_prot))
        chain_data.append(current_prot)






cols = [f.__name__ for f in features_to_compute]
# feature_avg_df = pandas.read_csv("feature_averages.csv")[cols]
# avg_dict = feature_avg_df.iloc[0,:].to_dict() 

# empty_row  = {k:0 for k in cols}
rows = []
# rows_raw = []
for prot in chain_data:
    if prot==[]:
        row  = avg_dict.copy()
        # rows_raw.append(empty_row.copy())
    else:
        row = calculate_features(prot)
        # rows_raw.append(row)
        for k in row.keys():
            try:
                a = row[k]/2
            except:
                row[k] = avg_dict[k]               
    rows.append(row)




# features_df_raw = pandas.DataFrame(rows_raw)

# features_df_raw.to_csv("features_temp_raw.csv")


features_df = pandas.DataFrame(rows)

# features_df.to_csv("features_temp.csv")
# features_df = features_df.fillna(feature_avg_df)



# features_df.to_csv(this_script_directory+"/"+taskid+"_chains.csv")

# with open(this_script_directory+"/"+taskid+"_inputdf.txt","w") as cmdfile:
#     cmdfile.writelines(str(features_df))
# out_df.to_csv(in_path[:-4]+"_features.txt",index = False, sep = "\t")

with open(this_script_directory +"/models/resolution.model","rb") as res_modelf,\
 		open(this_script_directory +"/models/rfree.model","rb") as rfree_modelf,\
 		open(this_script_directory +"/models/resolution.scaler","rb") as res_scalerf,\
        open(this_script_directory +"/models/rfree.scaler","rb") as rfree_scalerf,\
 		open(this_script_directory +"/conf.json","rb") as conff:
    res_model = load(res_modelf)
    rfree_model = load(rfree_modelf)
    res_scaler= load(res_scalerf)
    rfree_scaler = load(rfree_scalerf)
    all_conf= json.load(conff)





with open(this_script_directory +"/data/trainset_res_vals_sorted.txt") as train_resf , open(this_script_directory +"/data/trainset_rfree_vals_sorted.txt") as train_rfreef:
    train_res_vals_str = train_resf.read().rstrip("\n").split("\n")
    train_res_vals = [float(x) for x in train_res_vals_str]
    train_rfree_vals_str = train_rfreef.read().rstrip("\n").split("\n")
    train_rfree_vals = [float(x) for x in train_rfree_vals_str]

################ for resolution
conf= all_conf["resolution"]
min_val= train_res_vals[0]
max_val= train_res_vals[-1]
features = conf["features"]
norm_df = normalize_by_scaler(features_df,features,res_scaler)
pred_res = predict(res_model,norm_df,features,min_val,max_val)

################ for rfree
conf= all_conf["rfree"]
min_val= train_rfree_vals[0]
max_val= train_rfree_vals[-1]
features = conf["features"]
norm_df = normalize_by_scaler(features_df, features, rfree_scaler)
pred_rfree = predict(rfree_model,norm_df,features,min_val,max_val)



in_dir = os.path.dirname(in_path)
html_path = in_dir+ "/results.html"
csv_path = in_dir + "/results.csv"

# with open(this_script_directory+"/"+taskid+"_result.txt","w") as cmdfile:
#     cmdfile.writelines(html_path+"\n"+str(pids))

prepare_output(html_path, csv_path, pids, pred_res, pred_rfree,profile_exit_code,profile_stop_chain, train_res_vals, train_rfree_vals)
