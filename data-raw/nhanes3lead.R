## code to prepare `nhanes3lead` dataset goes here

library(tidyverse)

youth_cols <- list(
  SEQN = c(1,5), FAMILY_SEQN = c(6,10), STATUS = c(11,11), RACE_ETH = c(12,12),
  RACE = c(13,13), ETH = c(14,14), SEX = c(15,15), AGE = c(21,24), PIR = c(36,41),
  WEIGHTS = c(60,67),
  SMOKE_HOME = c(1292,1292), FOOD = c(1313,1313), EDU = c(1359,1360),
  SMOKE_PREG = c(1379,1379), BABY_NICU = c(1383,1383)
)
exam_cols <- list(
  SEQN = c(1,5), FAMILY_SEQN = c(6,10), STATUS = c(11,11),
  MATH = c(4433,4434), READING = c(4435,4436), BLOCK = c(4437,4438), DIGIT = c(4439,4440)
)
lab_cols <- list(
  SEQN = c(1,5), FAMILY_SEQN = c(6,10), STATUS = c(11,11), BLL = c(1423,1426)
)

read_ascii <- function(file, col_list) {
  read_fwf(file, fwf_positions(
    start = map_int(col_list, 1), end = map_int(col_list, 2),
    col_names = names(col_list)
  ), col_types = cols(.default = "c"))
}

raw_youth_df <- read_ascii("https://wwwn.cdc.gov/nchs/data/nhanes3/1a/youth.dat", youth_cols)
raw_exam_df <- read_ascii("https://wwwn.cdc.gov/nchs/data/nhanes3/1a/exam.dat", exam_cols)
raw_lab_df  <- read_ascii("https://wwwn.cdc.gov/nchs/data/nhanes3/1a/lab.dat",  lab_cols)

# 2. Merge data
raw_merged_df <- raw_youth_df |>
  full_join(select(raw_exam_df, -FAMILY_SEQN, -STATUS), by = "SEQN") |>
  full_join(select(raw_lab_df,  -FAMILY_SEQN, -STATUS), by = "SEQN")

# 3. Clean and recode data
raw_merged_df <- raw_merged_df |>
  mutate(
    STATUS = (recode(STATUS, "1"="No exam", "2"="MEC exam", "3"="Home exam")),
    RACE_ETH = (recode(RACE_ETH, "1"="White", "2"="Black", "3"="Hispanic", "4"="Other")),
    RACE = (recode(RACE, "1"="White", "2"="Black", "3"="Other", "8"="Mex_Am")),
    ETH = (recode(ETH, "1"="Mex Am", "2"="Other Hisp", "3"="Not Hisp")),
    SEX = (recode(SEX, "1"="Male", "2"="Female")),
    AGE = na_if(as.numeric(AGE), 9999) / 12,
    PIR = PIR |> as.numeric() |> na_if(888888),
    EDU = as.integer(EDU) |> na_if(88) |> na_if(99),
    SMOKE_HOME = recode(SMOKE_HOME, "1"="Yes", "2"="No", "8"="None"),
    FOOD = recode(FOOD, "1"="Good", "2"="Sometimes bad", "3"="Often bad", "8"="None"),
    SMOKE_PREG = recode(SMOKE_PREG, "1"="Yes", "2"="No", .default="None"),
    BABY_NICU = recode(BABY_NICU, "1" = "Yes", "2"="No", .default="None"),
    # Score columns: remove blanks/"88"/NA/"NaN"
    MATH = str_trim(MATH) |> na_if("88") |> na_if("") |> as.numeric(),
    READING = str_trim(READING) |> na_if("88") |> na_if("") |> as.numeric(),
    BLOCK = str_trim(BLOCK) |> na_if("88") |> na_if("") |> as.numeric(),
    DIGIT = str_trim(DIGIT) |> na_if("88") |> na_if("") |> as.numeric(),
    BLL = BLL |> na_if("8888") |> na_if("") |> as.numeric(),
    WEIGHTS = as.numeric(WEIGHTS)
  )

raw_merged_df <- raw_merged_df |>
  filter(
    !is.na(MATH), !is.na(READING), !is.na(BLOCK), !is.na(DIGIT),
    !is.na(BLL), !is.na(PIR),
    SMOKE_HOME != "None", FOOD != "None", SMOKE_PREG != "None", BABY_NICU != "None"
  )

# 4. Select and make dummies
nhanes3lead <- raw_merged_df |>
  select(RACE_ETH, SEX, AGE, PIR, SMOKE_HOME, FOOD, SMOKE_PREG, BABY_NICU,
         MATH, READING, BLOCK, DIGIT, BLL, WEIGHTS) |>
  filter(BLL <= 30) |> # only BLL <= 30
  mutate(across(where(is.character), as.factor)) |>
  rename(Math = MATH, Reading = READING, Block = BLOCK, Digit = DIGIT, Age = AGE,
         Race = RACE_ETH, Sex = SEX, Smoke_in_Home = SMOKE_HOME,
         Smoke_Pregnant = SMOKE_PREG, NICU = BABY_NICU, Enough_Food = FOOD,
         MEC_wt = WEIGHTS) |>
  mutate(Male = as.integer(Sex == "Male"),
         Smoke_in_Home = as.integer(Smoke_in_Home == "Yes"),
         Smoke_Pregnant = as.integer(Smoke_Pregnant == "Yes"),
         NICU = as.integer(NICU == "Yes"),
         Enough_Food = as.integer(Enough_Food == "Good"),
         logBLL = log(BLL)) |>
  select(-Sex) |>
  select(Age, Male, Race, PIR, Enough_Food, Smoke_in_Home, Smoke_Pregnant, NICU,
         logBLL, Math, Reading, Block, Digit, MEC_wt) |>
  as.data.frame()

usethis::use_data(nhanes3lead, overwrite = TRUE)
