pro mk_str

cat_in = {starnm:'', tag_in:'amp', tag_out:'amp60',templ_obnm:['n', 'n', 'n', 'n', 'n'], comment:''} 
cat = replicate(cat_in,100)

cat[0].starnm='166'     & cat[0].comment='Drop BY Dra'
                          cat[0].templ_obnm[0:2]='191008.'+['1074','1078', '1080'] 
cat[1].starnm='3651'    & cat[1].templ_obnm[0:1]='191008.'+['1084','1085']  
cat[2].starnm='4614'    & cat[2].comment='Wait Only 2018 obs'
                        ;  cat[2].templ_obnm[0:2]='.'+[]  
cat[3].starnm='4628'    & cat[3].templ_obnm[0:1]='191008.'+['1089','1090']  
cat[4].starnm='6582'    & cat[4].comment='Wait One night 2018'
                        ; cat[4].templ_obnm[0:2]='.'+[]  
cat[5].starnm='9407'    & cat[5].comment='Wait Need high SNR templ'
                        ; cat[5].templ_obnm[0:2]='.'+[]  
cat[6].starnm='10476'   & cat[6].templ_obnm[0:2]='191007.'+['1124','1125','1126']  
cat[7].starnm='10700'   & cat[7].templ_obnm[0:3]='191008.'+['1094','1095','1096','1097']  
cat[8].starnm='10780'   & cat[8].templ_obnm[0:1]='190821.'+['1117','1118']  
cat[9].starnm='16160'   & cat[9].comment='Wait Only 2 obs, 2018'
                        ; cat[9].templ_obnm[0:2]='.'+[]  
cat[10].starnm='17156'  & cat[10].comment='Wait Need high SNR templ'
                        ; cat[10].templ_obnm[0:2]='.'+[]  
cat[11].starnm='18803'  & cat[11].comment='Wait 1 obs'
                        ; cat[11].templ_obnm[0:2]='.'+[]  
cat[12].starnm='22049'  & cat[12].templ_obnm[0:4]='191008.'+['1125','1126','1127','1128','1129']  
cat[13].starnm='23249'  & cat[13].comment='Wait Only 2018 obs'
                        ; cat[13].templ_obnm[0:2]='.'+[]  
cat[14].starnm='25680'  & cat[14].templ_obnm[0:1]='191008.'+['1157','1158']  
cat[15].starnm='26965'  & cat[15].templ_obnm[0:2]='191008.'+['1133','1134','1135']  
cat[16].starnm='34411'  & cat[16].templ_obnm[0:2]='191008.'+['1168','1169','1170']  
cat[17].starnm='37394'  & cat[17].comment='Wait Only 2 nights'
                          cat[17].templ_obnm[0:1]='191008.'+['1174','1175']  
cat[18].starnm='50692'  & cat[18].templ_obnm[0:2]='191015.'+['1184','1185','1186']  
cat[19].starnm='62044'  & cat[19].comment='Drop rapid rotator' 
                          cat[19].templ_obnm[0:2]='190225.'+['1111','1112','1113']  
cat[20].starnm='72905'  & cat[20].comment='Drop rapid rotator'
                          cat[20].templ_obnm[0:2]='190210.'+['1108','1112','1116']  
cat[21].starnm='75732'  & cat[21].templ_obnm[0:3]='190225.'+['1079','1080','1084','1085']  
cat[22].starnm='80606'  & cat[22].templ_obnm[0:4]='190311.'+['1069','1070','1071','1072','1073']  
cat[23].starnm='81009'  & cat[23].comment='Wait One night 2018'
                        ; cat[23].templ_obnm[0:2]='.'+[]  
cat[24].starnm='86728'  & cat[24].templ_obnm[0:2]='190210.'+['1119','1120','1123']  
cat[25].starnm='88133'  & cat[25].comment='Wait One night 2018'
                        ; cat[25].templ_obnm[0:2]='.'+[]  
cat[26].starnm='88230'  & cat[26].comment='Drop SB2' 
                        ; cat[26].templ_obnm[0:1]='190419.'+['1106','1110']  
cat[27].starnm='89269'  & cat[27].comment='Wait Need high SNR template' 
                        ; cat[27].templ_obnm[0:2]='.'+[]  
cat[28].starnm='89744'  & cat[28].comment='Known eccentric Jupiter; rapid rotator'    
                          cat[28].templ_obnm[0:1]='190315.'+['1081','1082']  
cat[29].starnm='95128'  & cat[29].comment='47 UMa - known planets'
                          cat[29].templ_obnm[0:2]='190210.'+['1126','1127','1129']  
cat[30].starnm='95735'  & cat[30].templ_obnm[0:1]='190419.'+['1119','1123']  
cat[31].starnm='99491'  & cat[31].templ_obnm[0:1]='190503.'+['1077','1078']  
cat[32].starnm='99492'  & cat[32].templ_obnm[0:2]='190518.'+['1113','1117','1118']  
cat[33].starnm='99995'  & cat[33].templ_obnm[0:2]='190210.'+['1132','1135','1138']  
cat[34].starnm='101177' & cat[34].comment='Wait few obs in 2018 only' 
                        ; cat[34].templ_obnm[0:2]='.'+[]  
cat[35].starnm='101501' & cat[35].templ_obnm[0:2]='190210.'+['1153','1154','1155']  
cat[36].starnm='103095' & cat[36].comment='Wait Need high SNR template obs' 
                        ; cat[36].templ_obnm[0:2]='.'+[]  
cat[37].starnm='105631' & cat[37].templ_obnm[0:1]='190503.'+['1092','1093']  
cat[38].starnm='109358' & cat[38].templ_obnm[0:2]='190210.'+['1158','1159','1160']  
cat[39].starnm='114783' & cat[39].templ_obnm[0:2]='190503.'+['1098','1099','1100']  
cat[40].starnm='115617' & cat[40].templ_obnm[0:2]='190210.'+['1163','1164','1165']  
cat[41].starnm='117043' & cat[41].templ_obnm[0:1]='190503.'+['1104','1105']  
cat[42].starnm='120066' & cat[42].comment='Drop - SB?' 
                          cat[42].templ_obnm[0:1]='190503.'+['1109','1110']  
cat[43].starnm='120136' & cat[43].comment='Drop rapid rotator' 
                          cat[43].templ_obnm[0:3]='190315.'+['1092','1093','1094','1095']  
cat[44].starnm='122064' & cat[44].templ_obnm[0:1]='190503.'+['1114','1115']  
cat[45].starnm='126053' & cat[45].comment='Wait Need high SNR templ' 
                        ;  cat[45].templ_obnm[0:2]='.'+[]  
cat[46].starnm='127334' & cat[46].templ_obnm[0:1]='190210.'+['1169','1170']  
cat[47].starnm='131511' & cat[47].comment='Wait Need high SNR templ' 
                        ; cat[47].templ_obnm[0:2]='.'+[]  
cat[48].starnm='135599' & cat[48].comment='Wait Need high SNR templ' 
                        ;cat[48].templ_obnm[0:2]='.'+[]  
cat[49].starnm='136923' & cat[49].comment='Wait Only 2018 obs' 
                        ;  cat[49].templ_obnm[0:2]='.'+[]  
cat[50].starnm='140538' & cat[50].comment='Wait Need high SNR templ' 
                        ; cat[50].templ_obnm[0:2]='.'+[]  
cat[51].starnm='141004' & cat[51].comment='Wait Only 2018 observations' 
                        ;  cat[51].templ_obnm[0:2]='.'+[]  
cat[52].starnm='143761' & cat[52].comment='Wait Only 3 obs' 
                        ; cat[52].templ_obnm[0:2]='.'+[]  
cat[53].starnm='144579' & cat[53].templ_obnm[0:1]='190427.'+['1085','1086']  
cat[54].starnm='145148' & cat[54].comment='Wait Only 2018 obs' 
                        ; cat[54].templ_obnm[0:2]='.'+[]  
cat[55].starnm='145675' & cat[55].templ_obnm[0:1]='190505.'+['1107','1108']  
cat[56].starnm='146233' & cat[56].comment='Wait Need high snr template' 
                        ; cat[56].templ_obnm[0:2]='.'+[]  
cat[57].starnm='149661' & cat[57].comment='Wait Need high SNR templ' 
                        ;  cat[57].templ_obnm[0:2]='.'+[]  
cat[58].starnm='152391' & cat[58].comment='Wait Only 2018 obs' 
                        ; cat[58].templ_obnm[0:2]='.'+[]  
cat[59].starnm='154345' & cat[59].comment='Wait Only 2 obs' 
                        ; cat[59].templ_obnm[0:2]='.'+[]  
cat[60].starnm='157214' & cat[60].templ_obnm[0:1]='190828.'+['1222','1223']  
cat[61].starnm='157347' & cat[61].templ_obnm[0:1]='190525.'+['1088','1089']  
cat[62].starnm='158259' & cat[62].comment='Wait Only 2 obs'
                        ; cat[62].templ_obnm[0:2]='.'+[]  
cat[63].starnm='158633' & cat[63].templ_obnm[0:1]='190817.'+['1073','1074']  
cat[64].starnm='159062' & cat[64].comment='Wait Only 2 obs' 
                        ; cat[64].templ_obnm[0:2]='.'+[]  
cat[65].starnm='159222' & cat[65].comment='Wait Only 2018 obs' 
                        ; cat[65].templ_obnm[0:2]='.'+[]  
cat[66].starnm='161797' & cat[66].templ_obnm[0:4]='190923.'+['1080','1081','1082','1083','1084']  
cat[67].starnm='164922' & cat[67].comment='Wait Need high SNR templ'
                        ; cat[67].templ_obnm[0:2]='.'+[]  
cat[68].starnm='165341' & cat[68].comment='Wait Need high SNR templ'
                        ; cat[68].templ_obnm[0:2]='.'+[]  
cat[69].starnm='166620' & cat[69].templ_obnm[0:1]='190817.'+['1078','1082']  
cat[70].starnm='168009' & cat[70].templ_obnm[0:1]='190828.'+['1239','1241']  
cat[71].starnm='182488' & cat[71].templ_obnm[0:1]='190831.'+['1108','1112']  
cat[72].starnm='183123' & cat[72].comment='Wait Only 2018 obs' 
                        ; cat[72].templ_obnm[0:2]='.'+[]  
cat[73].starnm='183756' & cat[73].comment='Wait Only 2018 obs' 
                        ; cat[73].templ_obnm[0:2]='.'+[]  
cat[74].starnm='185144' & cat[74].templ_obnm[0:2]='191007.'+['1070','1071','1072']  
cat[75].starnm='185603' & cat[75].comment='Wait Only 2018 observations'
                        ; cat[75].templ_obnm[0:2]='.'+[]  
cat[76].starnm='186408' & cat[76].templ_obnm[0:1]='190817.'+['1094','1095']  
cat[77].starnm='186427' & cat[77].templ_obnm[0:1]='190816.'+['1084','1085']  
cat[78].starnm='189733' & cat[78].templ_obnm[0:3]='190524.'+['1090','1091','1095','1096']  
cat[79].starnm='190404' & cat[79].comment='Wait Only 2018 obs' 
                        ; cat[79].templ_obnm[0:2]='.'+[]  
cat[80].starnm='190406' & cat[80].templ_obnm[0:2]='190824.'+['1293','1294','1295']  
cat[81].starnm='191785' & cat[81].comment='Wait Only 2018 obs'
                        ; cat[81].templ_obnm[0:2]='.'+[]  
cat[82].starnm='193664' & cat[82].comment='Wait Only 2018 obs'
                        ; cat[82].templ_obnm[0:2]='.'+[]  
cat[83].starnm='195689' & cat[83].comment='Wait Only 1 night'
                        ; cat[83].templ_obnm[0:2]='.'+[]  
cat[84].starnm='197076' & cat[84].templ_obnm[0:1]='190909.'+['1086','1091']  
cat[85].starnm='199960' & cat[85].comment='Wait One obs'
                        ; cat[85].templ_obnm[0:2]='.'+[]  
cat[86].starnm='201091' & cat[86].templ_obnm[0:2]='191007.'+['1076','1077','1078']  
cat[87].starnm='201092' & cat[87].comment='Wait One obs'
                        ; cat[87].templ_obnm[0:2]='.'+[]  
cat[88].starnm='210277' & cat[88].templ_obnm[0:1]='190923.'+['1121','1122']  
cat[89].starnm='217014' & cat[89].templ_obnm[0:1]='191008.'+['1070','1073']  
cat[90].starnm='218868' & cat[90].templ_obnm[0:1]='190819.'+['1105','1106']  
cat[91].starnm='219134' & cat[91].templ_obnm[0:2]='191007.'+['1087','1088','1093']  
cat[92].starnm='221354' & cat[92].templ_obnm[0:1]='190923.'+['1142','1146']  
cat[93].starnm=''  ;& cat[93].templ_obnm[0:2]='.'+[]  
cat[94].starnm=''  ;& cat[94].templ_obnm[0:2]='.'+[]  
cat[95].starnm=''  ;& cat[95].templ_obnm[0:2]='.'+[]  
cat[96].starnm=''  ;& cat[96].templ_obnm[0:2]='.'+[]  
cat[97].starnm=''  ;& cat[97].templ_obnm[0:2]='.'+[]  
cat[98].starnm=''  ;& cat[98].templ_obnm[0:2]='.'+[]  
cat[99].starnm=''  ;& cat[99].templ_obnm[0:2]='.'+[]  

save, cat, f='cat_drive.dat' 

stop
end
