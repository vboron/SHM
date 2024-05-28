#!/usr/bin/env python3

import pandas as pd
import runAGL

runAGL.shm_graphing(pd.read_csv('free_mutations.csv'), pd.read_csv('complex_mutations.csv'), pd.read_csv('complex_free_mutations.csv'))