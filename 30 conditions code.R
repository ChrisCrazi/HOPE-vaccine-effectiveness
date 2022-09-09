alcohol_misuse <- all_dx %>%
  filter(str_detect(code, "^265.2|^291.1|^291.2|^291.3|^291.5|^291.6|^291.7|
                           ^291.8|^291.9|^303.0|^303.9|^305.0|^357.5|^425.5|
                           ^535.3|^571.0|^571.1|^571.2|^571.3|^980|^V11.3|^P15")) %>% 
  mutate(dx = "Alcohol misuse")

asthma <- all_dx %>% 
  filter(str_detect(code, "^493|R96")) %>% 
  mutate(dx = "Asthma")

cancer_lymphoma <- all_dx %>% 
  filter(str_detect(code, "^200|^201|^202|^203.0|^238.6|^B72")) %>% 
  mutate(dx = "Cancer, lymphoma")

cancer_metastatic <- all_dx %>% 
  filter(str_detect(code, "^196|^197|^198|^199|^B74|^D74|^D76|^D77|^L71|^N74|^S77|^T71|^U75|^U76|
                           ^U77|^W72|^X77|^Y78")) %>% 
  mutate(dx = "Cancer, metastatic")

cancer_non_metastatic <- all_dx %>% 
  filter(str_detect(code, "^153|^154|^162|^163|^174|^180|185^|^230.3|^230.4|
                             ^230.5|^230.6|^231.2|^233.0|^233.1|^233.4|^D75|^R84|^X75|
                             ^X76|^Y77")) %>% 
  mutate(dx = "Cancer, non-metastatic")

chronic_pain <- all_dx %>% 
  filter(str_detect(code, "^307.80|^307.89|^338.0|^338.2|^338.4|^719.41|
                             ^719.45|^719.46|^719.47|^719.49|^720.0|^720.2|
                             ^720.9|^721.0|^721.1|^721.2|^721.3|^721.4|^721.6|^721.8|^721.9|
                             ^722|^723.0|^723.1|^723.[3456789]|^724.0|^724.1|
                             ^724.2|^724.3|^724.4|^724.5|^724.6|^724.70|^724.79|
                             ^724.8|^724.9|^729.0|^729.1|^729.2|^729.4|^729.5|^A01")) %>% 
  mutate(dx = "Chronic pain")

chronic_pulmonary_disease <- all_dx %>% 
  filter(str_detect(code, "^416.8|^416.9|^490|^491|^492|^494|^495|^496|^497|^498|
                             ^499|^500|^501|^502|^503|^504|^505|^506.4|^508.1|^508.8|^R95")) %>% 
  mutate(dx = "Chronic pulmonary disease")

chronic_viral_hepatitis_b <- all_dx %>% 
  filter(str_detect(code, "^70.2|^70.3|^D72")) %>% 
  mutate(dx = "Chronic viral hepatitis B")

cirrhosis <- all_dx %>% 
  filter(str_detect(code, "^571.2|^571.5|^571.6|^456.0|^456.1|^456.20|^456.21|^567.0|^567.2|
                             ^567.21|^567.29|^567.8|^567.9|^572.2|^572.4|^789.5")) %>% 
  mutate(dx = "Cirrhosis")

dementia <- all_dx %>% 
  filter(str_detect(code, "^290|^294.1|^331.2|^P70")) %>% 
  mutate(dx = "Dementia")

depression <- all_dx %>% 
  filter(str_detect(code, "^296.2|^296.3|^296.5|^300.4|^309|^311|^P76")) %>% 
  mutate(dx = "Depression")

diabetes <- all_dx %>% 
  filter(str_detect(code, "^250|^T89|^T90")) %>% 
  mutate(dx = "Diabetes")

hypertension <- all_dx %>% 
  filter(str_detect(code, "^401|^402|^403|^404|^405|^K86|^K87")) %>% 
  mutate(dx = "Hypertension")

hypothyroidism <- all_dx %>% 
  filter(str_detect(code, "^240.9|^243|^244|^246.1|^246.8|^T86")) %>% 
  mutate(dx = "Hypothyroidism")

inflammatory_bowel_disease <- all_dx %>% 
  filter(str_detect(code, "^555|^556")) %>% 
  mutate(dx = "Inflammatory bowel disease")

irritable_bowel_syndrome <- all_dx %>% 
  filter(str_detect(code, "^564.1|^D93")) %>% 
  mutate(dx = "Irritable bowel syndrome")

parkinson <- all_dx %>% 
  filter(str_detect(code, "^332|^N87")) %>% 
  mutate(dx = "Parkinson's disease")

peptic_ulcer_disease <- all_dx %>% 
  filter(str_detect(code, "^531.7|^531.9|^532.7|^532.9|^533.7|^533.9|^534.7|^534.9|^D86")) %>% 
  mutate(dx = "Peptic ulcer disease")

peripheral_vascular_disease <- all_dx %>% 
  filter(str_detect(code, "^440.2|^K92")) %>% 
  mutate(dx = "Peripheral vascular disease")

psoriasis <- all_dx %>% 
  filter(str_detect(code, "^696.1|^S91")) %>% 
  mutate(dx = "Psoriasis")

rheumatoid_arthritis <- all_dx %>% 
  filter(str_detect(code, "^446.5|^710.0|^710.1|^710.2|^710.3|^710.4|^714.0|^714.1|
                           ^714.2|^714.8|^725|^L88")) %>% 
  mutate(dx = "Rheumatoid arthritis")

schizophrenia <- all_dx %>% 
  filter(str_detect(code, "^295|^P72")) %>% 
  mutate(dx = "Schizophrenia")

severe_constipation <- all_dx %>% 
  filter(str_detect(code, "^560.1|^560.30|^560.39|^560.9|^564.0|^569.83|^569.89|^D12")) %>% 
  mutate(dx = "Severe constipation")

#-----------------------------------------------------

atrial_fibrillation <- all_dx %>%
  filter(str_detect(code, "^427.3|^K78")) %>% 
  mutate(dx = "Atrial fibrillation")

chronic_heart_failure <- all_dx %>% 
  filter(str_detect(code, "^398.91|^402.01|^402.11|^402.91|^404.01|^404.03|
                             ^404.11|^404.13|^404.91|^404.93|^425.4|^425.5|
                             ^425.6|^425.7|^425.8|^425.9|^428|^K77")) %>% 
  mutate(dx = "Chronic heart failure")


chronic_kidney_disease <- all_dx %>% 
  filter(str_detect(code, "^583|^584|^585|^586|^592|^593.9|^U14")) %>% 
  mutate(dx = "Chronic kidney disease")

epilepsy <- all_dx %>% 
  filter(str_detect(code, "^345|^N88")) %>% 
  mutate(dx = "Epilepsy")

multiple_sclerosis <- all_dx %>% 
  filter(str_detect(code, "^323|^340|^341.0|^341.9|^377.3|^N86")) %>% 
  mutate(dx = "Multiple sclerosis")

myocardial_infarction <- all_dx %>% 
  filter(str_detect(code, "^410|^K75")) %>% 
  mutate(dx = "Myocardial infarction")

stroke_or_tia <- all_dx %>% 
  filter(str_detect(code, "^326.3|^430|^431|^433.[0123456789]1|^434.[0123456789]1|^435|^436|^K90")) %>% 
  mutate(dx = "Stroke or TIA")

all_morbidity_dx <- rbind(alcohol_misuse, asthma, cancer_lymphoma, cancer_metastatic, cancer_non_metastatic,
                          chronic_pain, chronic_pulmonary_disease, chronic_viral_hepatitis_b,
                          cirrhosis, dementia, depression, diabetes, hypertension, hypothyroidism, 
                          inflammatory_bowel_disease, irritable_bowel_syndrome, parkinson,
                          peptic_ulcer_disease, peripheral_vascular_disease, psoriasis, rheumatoid_arthritis,
                          schizophrenia, severe_constipation,
                          atrial_fibrillation,
                          chronic_heart_failure,
                          chronic_kidney_disease,
                          epilepsy,
                          multiple_sclerosis,
                          myocardial_infarction,
                          stroke_or_tia)

mm_pssn <- all_morbidity_dx %>% 
  pull(patient_pssn) %>% 
  unique()
