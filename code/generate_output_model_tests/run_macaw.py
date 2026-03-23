import cobra
from macaw.main import run_all_tests, dead_end_test, duplicate_test, dilution_test, loop_test

model = cobra.io.read_sbml_model('model/iDG1237.xml')
#open up all of the exchange reactions:
for reac in model.reactions:
    if reac.id.startswith('EX_'):
        reac.bounds = (-1000, 1000)

# (test_results, edge_list) = run_all_tests(model)

# test_result.to_csv(out_file, index=False)

(dead_end_results, dead_end_edges) = dead_end_test(model)
(duplicate_results, duplicate_edges) = duplicate_test(model)
output = dead_end_results.merge(duplicate_results)
out_file = "data/Test_results/macaw_report_iDG1237.csv"
output.to_csv(out_file, index = False)

# NB: there is potential to run more tests, maybe try to run them all --> so far, it does not work.

# (dilution_results, dilution_edges) = dilution_test(model, dead_end_results=dead_end_results)
# # output = output.merge(dilution_results)
#
# # #collect all metabolite ID's for diphosphate/phosphate
# # (diphosphate_results, diphosphate_edges) = diphosphate_test(model, [], [])
#
# # (loop_results, loop_edges) = loop_test(model)
# # output = output.merge(loop_results)
#
# print(dilution_results.head())
