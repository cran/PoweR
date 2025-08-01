/*
Note: stat84 is available and can be used to incorporate new tests
 */
static const R_CMethodDef CEntries[] = {
  {"calcfx",        (DL_FUNC) &calcfx,       5},
  {"matrixpval",   (DL_FUNC) &matrixpval,   13},
  {"matrixpvalMC", (DL_FUNC) &matrixpvalMC, 17},
  {"compquantc",    (DL_FUNC) &compquantc,    17},
  {"powcompeasy",  (DL_FUNC) &powcompeasy,  13},
  {"powcompfast",  (DL_FUNC) &powcompfast,  29},
  {"statcomputeC",  (DL_FUNC) &statcomputeC,  17},
  // HERE WE PUT THE lawj NAMES
  {"law1",   (DL_FUNC) &law1,  7},
  {"law2",   (DL_FUNC) &law2,  7},
  {"law3",   (DL_FUNC) &law3,  7},
  {"law4",   (DL_FUNC) &law4,  7},
  {"law5",   (DL_FUNC) &law5,  7},
  {"law6",   (DL_FUNC) &law6,  7},
  {"law7",   (DL_FUNC) &law7,  7},
  {"law8",   (DL_FUNC) &law8,  7},
  {"law9",   (DL_FUNC) &law9,  7},
  {"law10",  (DL_FUNC) &law10, 7},
  {"law11",  (DL_FUNC) &law11, 7},
  {"law12",  (DL_FUNC) &law12, 7},
  {"law13",  (DL_FUNC) &law13, 7},
  {"law14",  (DL_FUNC) &law14, 7},
  {"law15",  (DL_FUNC) &law15, 7},
  {"law16",  (DL_FUNC) &law16, 7},
  {"law17",  (DL_FUNC) &law17, 7},
  {"law18",  (DL_FUNC) &law18, 7},
  {"law19",  (DL_FUNC) &law19, 7},
  {"law20",  (DL_FUNC) &law20, 7},
  {"law21",  (DL_FUNC) &law21, 7},
  {"law22",  (DL_FUNC) &law22, 7},
  {"law23",  (DL_FUNC) &law23, 7},
  {"law24",  (DL_FUNC) &law24, 7},
  {"law25",  (DL_FUNC) &law25, 7},
  {"law26",  (DL_FUNC) &law26, 7},
  {"law27",  (DL_FUNC) &law27, 7},
  {"law28",  (DL_FUNC) &law28, 7},
  {"law29",  (DL_FUNC) &law29, 7},
  {"law30",  (DL_FUNC) &law30, 7},
  {"law31",  (DL_FUNC) &law31, 7},
  {"law32",  (DL_FUNC) &law32, 7},
  {"law33",  (DL_FUNC) &law33, 7},
  {"law34",  (DL_FUNC) &law34, 7},
  {"law35",  (DL_FUNC) &law35, 7},
  {"law36",  (DL_FUNC) &law36, 7},
  {"law37",  (DL_FUNC) &law37, 7},
  {"law38",  (DL_FUNC) &law38, 7},
  {"law39",  (DL_FUNC) &law39, 7},
  {"law40",  (DL_FUNC) &law40, 7},
  {"law41",  (DL_FUNC) &law41, 7},
  {"law42",  (DL_FUNC) &law42, 7},
  {"law43",  (DL_FUNC) &law43, 7},
  // HERE WE PUT THE statj NAMES
  {"stat1",  (DL_FUNC) &stat1,  16},
  {"stat2",  (DL_FUNC) &stat2,  16},
  {"stat3",  (DL_FUNC) &stat3,  16},
  {"stat4",  (DL_FUNC) &stat4,  16},
  {"stat5",  (DL_FUNC) &stat5,  16},
  {"stat6",  (DL_FUNC) &stat6,  16},
  {"stat7",  (DL_FUNC) &stat7,  16},
  {"stat8",  (DL_FUNC) &stat8,  16},
  {"stat9",  (DL_FUNC) &stat9,  16},
  {"stat10", (DL_FUNC) &stat10, 16},
  {"stat11", (DL_FUNC) &stat11, 16},
  {"stat12", (DL_FUNC) &stat12, 16},
  {"stat13", (DL_FUNC) &stat13, 16},
  {"stat14", (DL_FUNC) &stat14, 16},
  {"stat15", (DL_FUNC) &stat15, 16},
  {"stat16", (DL_FUNC) &stat16, 16},
  {"stat17", (DL_FUNC) &stat17, 16},
  {"stat18", (DL_FUNC) &stat18, 16},
  {"stat19", (DL_FUNC) &stat19, 16},
  {"stat20", (DL_FUNC) &stat20, 16},
  {"stat21", (DL_FUNC) &stat21, 16},
  {"stat22", (DL_FUNC) &stat22, 16},
  {"stat23", (DL_FUNC) &stat23, 16},
  {"stat24", (DL_FUNC) &stat24, 16},
  {"stat25", (DL_FUNC) &stat25, 16},
  {"stat26", (DL_FUNC) &stat26, 16},
  {"stat27", (DL_FUNC) &stat27, 16},
  {"stat28", (DL_FUNC) &stat28, 16},
  {"stat29", (DL_FUNC) &stat29, 16},
  {"stat30", (DL_FUNC) &stat30, 16},
  {"stat31", (DL_FUNC) &stat31, 16},
  {"stat32", (DL_FUNC) &stat32, 16},
  {"stat33", (DL_FUNC) &stat33, 16},
  {"stat34", (DL_FUNC) &stat34, 16},
  {"stat35", (DL_FUNC) &stat35, 16},
  {"stat36", (DL_FUNC) &stat36, 16},
  {"stat37", (DL_FUNC) &stat37, 16},
  {"stat38", (DL_FUNC) &stat38, 16},
  {"stat39", (DL_FUNC) &stat39, 16},
  {"stat40", (DL_FUNC) &stat40, 16},
  {"stat41", (DL_FUNC) &stat41, 16},
  {"stat42", (DL_FUNC) &stat42, 16},
  {"stat43", (DL_FUNC) &stat43, 16},
  {"stat44", (DL_FUNC) &stat44, 16},
  {"stat45", (DL_FUNC) &stat45, 16},
  {"stat46", (DL_FUNC) &stat46, 16},
  {"stat47", (DL_FUNC) &stat47, 16},
  {"stat48", (DL_FUNC) &stat48, 16},
  {"stat49", (DL_FUNC) &stat49, 16},
  {"stat50", (DL_FUNC) &stat50, 16},
  {"stat51", (DL_FUNC) &stat51, 16},
  {"stat52", (DL_FUNC) &stat52, 16},
  {"stat53", (DL_FUNC) &stat53, 16},
  {"stat54", (DL_FUNC) &stat54, 16},
  {"stat55", (DL_FUNC) &stat55, 16},
  {"stat56", (DL_FUNC) &stat56, 16},
  {"stat57", (DL_FUNC) &stat57, 16},
  {"stat58", (DL_FUNC) &stat58, 16},
  {"stat59", (DL_FUNC) &stat59, 16},
  {"stat60", (DL_FUNC) &stat60, 16},
  {"stat61", (DL_FUNC) &stat61, 16},
  {"stat62", (DL_FUNC) &stat62, 16},
  {"stat63", (DL_FUNC) &stat63, 16},
  {"stat64", (DL_FUNC) &stat64, 16},
  {"stat65", (DL_FUNC) &stat65, 16},
  {"stat66", (DL_FUNC) &stat66, 16},
  {"stat67", (DL_FUNC) &stat67, 16},
  {"stat68", (DL_FUNC) &stat68, 16},
  {"stat69", (DL_FUNC) &stat69, 16},
  {"stat70", (DL_FUNC) &stat70, 16},
  {"stat71", (DL_FUNC) &stat71, 16},
  {"stat72", (DL_FUNC) &stat72, 16},
  {"stat73", (DL_FUNC) &stat73, 16},
  {"stat74", (DL_FUNC) &stat74, 16},
  {"stat75", (DL_FUNC) &stat75, 16},
  {"stat76", (DL_FUNC) &stat76, 16},
  {"stat77", (DL_FUNC) &stat77, 16},
  {"stat78", (DL_FUNC) &stat78, 16},
  {"stat79", (DL_FUNC) &stat79, 16},
  {"stat80", (DL_FUNC) &stat80, 16},
  {"stat81", (DL_FUNC) &stat81, 16},
  {"stat82", (DL_FUNC) &stat82, 16},
  {"stat83", (DL_FUNC) &stat83, 16},
  {"stat84", (DL_FUNC) &stat84, 16},
  {"stat85", (DL_FUNC) &stat85, 16},
  {"stat86", (DL_FUNC) &stat86, 16},
  {"stat87", (DL_FUNC) &stat87, 16},
  {"stat88", (DL_FUNC) &stat88, 16},
  {"stat89", (DL_FUNC) &stat89, 16},
  {"stat90", (DL_FUNC) &stat90, 16},
  {"stat91", (DL_FUNC) &stat91, 16},
  {"stat92", (DL_FUNC) &stat92, 16},
  {"stat93", (DL_FUNC) &stat93, 16},
  {"stat94", (DL_FUNC) &stat94, 16},
  {"stat95", (DL_FUNC) &stat95, 16},
  {"stat96", (DL_FUNC) &stat96, 16},
  {"stat97", (DL_FUNC) &stat97, 16},
  {"stat98", (DL_FUNC) &stat98, 16},
  {"stat99", (DL_FUNC) &stat99, 16},
  {"stat100", (DL_FUNC) &stat100, 16},
  {"stat101", (DL_FUNC) &stat101, 16},
  {"stat102", (DL_FUNC) &stat102, 16},
  {"stat103", (DL_FUNC) &stat103, 16},
  {"stat104", (DL_FUNC) &stat104, 16},
  {"stat105", (DL_FUNC) &stat105, 16},
  {"stat106", (DL_FUNC) &stat106, 16},
  {"stat107", (DL_FUNC) &stat107, 16},
  {"stat108", (DL_FUNC) &stat108, 16},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"compquantRcpp",    (DL_FUNC) &compquantRcpp,    19},
    {"gensampleRcpp",    (DL_FUNC) &gensampleRcpp,     7},
    {"matrixpvalMCRcpp", (DL_FUNC) &matrixpvalMCRcpp, 20},
    {"matrixpvalRcpp",   (DL_FUNC) &matrixpvalRcpp,   15},
    {"powcompeasyRcpp",  (DL_FUNC) &powcompeasyRcpp,  16},
    {"powcompfastRcpp",  (DL_FUNC) &powcompfastRcpp,  32},
    {"statcomputeRcpp",  (DL_FUNC) &statcomputeRcpp,   3},
    {NULL, NULL, 0}
};

void R_init_PoweR(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
