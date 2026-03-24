import cobra
import json
import csv
import itertools
import pandas as pd

def json_load(path):
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: str ~ The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary


def get_enzyme_constraint_model(json_model_file):
    """using enzyme concentration constraint
    json model to create a COBRApy model.

    Arguments
    ----------
    * json_model_file: json Model file.

    :return: Construct an enzyme-constrained model.
    """

    dictionary_model = json_load(json_model_file)
    model = cobra.io.json.load_json_model(json_model_file)

    coefficients = dict()
    for rxn in model.reactions:
        for eachr in dictionary_model['reactions']:
            if rxn.id == eachr['id']:
                if eachr['kcat_MW']:
                    coefficients[rxn.forward_variable] = 1 / float(eachr['kcat_MW'])
                break

    lowerbound = dictionary_model['enzyme_constraint']['lowerbound']
    upperbound = dictionary_model['enzyme_constraint']['upperbound']
    constraint = model.problem.Constraint(0, lb=lowerbound, ub=upperbound, name='Proteome_constraint')
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)
    return model


def get_enzyme_constraint_sep_protein_pools(proteome_file, json_model_file, proteome_dict, verbose=False):
    """

    :param proteome_file:
    :param json_model_file:
    :param proteome_dict: keys should match the: 'STEPdb Sub-cellular Location (Letter Code)' column of the proteome file,
        should always have an 'other' category.
    :return:
    """
    #NB first try to adapt this code for a proteome file that is based on reaction-ids!!!
    link_reac_to_cat_dict = {}
    coefficients_dict = {}

    # for category in proteome_dict:
    #     reaction_dicts[category] = []
    category_list = list(proteome_dict.keys())
    for cat in proteome_dict.keys():
        coefficients_dict[cat] = {}
    with open(proteome_file, 'r') as proteome:
        proteome_reader = csv.DictReader(proteome, delimiter='\t')
        for row in proteome_reader:
            has_category = False
            row_cats = row['Pool'].split(',')
            reac_id = row['Reaction']
            link_reac_to_cat_dict[reac_id] = []
            for cat in row_cats:
                if cat in category_list:
                    has_category = True
                    link_reac_to_cat_dict[reac_id].append(cat)
            if not has_category:
                link_reac_to_cat_dict[reac_id].append('other')


    dictionary_model = json_load(json_model_file)
    model = cobra.io.json.load_json_model(json_model_file)
    for eachr in dictionary_model['reactions']:
        if eachr['kcat_MW']:
            reac_id = eachr['id']
            reaction_object = model.reactions.get_by_id(reac_id)
            try:
                cat_list = link_reac_to_cat_dict[reac_id]
                if verbose:
                    print(reac_id, cat_list)
            except KeyError:
                cat_list = ["other"]
            for cat in cat_list:
                coefficients_dict[cat][reaction_object.forward_variable] = 1 / float(eachr['kcat_MW'])

    for cat in coefficients_dict:
        coeff_dict = coefficients_dict[cat]
        if verbose:
            print(cat, len(coeff_dict))
        constraint_name = f"proteome_constraint_{cat}"
        lower_bound = 0
        upper_bound = proteome_dict[cat]
        new_constraint = model.problem.Constraint(0, lb=lower_bound, ub=upper_bound, name=constraint_name)
        model.add_cons_vars(new_constraint)
        model.solver.update()
        new_constraint.set_linear_coefficients(coefficients=coeff_dict)
    return model


def update_ub_proteome(model, new_ub_value):
    # take all of the (updated) coefficients and make a dict)
    const = model.constraints.Proteome_constraint
    const.ub = new_ub_value
    return model


def try_basic_enz_model_functions(model):
    enz_model = get_enzyme_constraint_model(model)
    #enz_model.objective = 'ATPM'

    solution = enz_model.optimize()
    print(solution)

    #enz_model.reactions.EX_pi_e.lower_bound = -0.14
    enz_model.reactions.EX_pi_e.lower_bound = -0.145
    #enz_model.reactions.EX_pi_e.lower_bound = -0.18
    enz_model.reactions.EX_nh4_e.bounds = 0, 1000
    #enz_model.reactions.EX_nh4_e_reverse.bounds = 0, 1000
    enz_model.reactions.SOYPEPT.bounds = 0, 0
    enz_model.reactions.Growth.lower_bound = 0.1

    #update_ub_proteome(enz_model, 0.227)

    update_ub_proteome(enz_model, 0.10)
    #update_ub_proteome(enz_model, 0.01)
    #
    # print(enz_model.reactions.PDH_num1.kcat_MW)
    #
    solution = enz_model.optimize()
    print(solution)
    for reac in enz_model.reactions:
        if reac.id.startswith('EX_'):
            flux = solution.fluxes[reac.id]
            if flux != 0:
                print(reac.id, flux)
    #print(enz_model.metabolites.nh4_c.summary(solution=solution))
    print(enz_model.metabolites.lys__L_c.summary(solution=solution))

    # enz_model.reactions.EX_gclavam_e.bounds = 0,0
    # solution = enz_model.optimize()
    # print(solution)
    # for reac in enz_model.reactions:
    #     if reac.id.startswith('EX_'):
    #         flux = solution.fluxes[reac.id]
    #         if flux != 0:
    #             print(reac.id, flux)

    return None


def quick_try_non_constrained_model(model_file):
    model = cobra.io.read_sbml_model(model_file)
    model.reactions.EX_pi_e.lower_bound = -0.14

    sol = model.optimize()
    print(sol)
    for reac in model.reactions:
        if reac.id.startswith('EX_'):
            flux = sol.fluxes[reac.id]
            if flux != 0:
                print(reac.id, flux)
    return None


def sample_model_with_pfba(model, settings):
    # copy the model so that only the copied model is changed by the pFBA constraints
    model_copy = model.copy()
    objective_reac = model_copy.reactions.get_by_id(settings['objective'])
    add_pFBA_constraint(model_copy, objective_reac, settings['objective_fraction'], settings['flux_fraction'])

    # fix the biomass reaction in the model before sampling!
    solution = model_copy.optimize()
    if 'fraction_of_optimum' in settings:
        frac_optimum = settings['fraction_of_optimum']
    else:
        frac_optimum = 0.95
    objective_reac.lower_bound = frac_optimum * solution.objective_value

    if 'threads' in settings:
        threads = settings['threads']
    else:
        threads = None

    # sampling part
    set_sampler = cobra.sampling.OptGPSampler(model=model_copy, processes=threads, thinning=1000, seed=123)
    sample_frame = set_sampler.sample(settings['n_samples'])
    return sample_frame


def add_pFBA_constraint(model, objective, objective_fraction=1, flux_fraction=1):
    """Add a pFBA flux minimality constraint to the model. Updates the model in-place.
    Function adapted from: https://gitlab.com/wurssb/Modelling/sampling-tools, subject to MIT license

    :param model: The cobra model
    :type model: cobra.Model
    :param objective: The objective that should be maximized,
    in order to determine the objective bound for the flux sum minimization.
    :type objective: cobra.Reaction
    :param objective_fraction: pFBA minimum objective fraction, should be 1<= and defaults to 1.
    :type objective_fraction: float, optional
    :param flux_fraction: pFBA maximum flux fraction, should be >=1 and defaults to 1.
    :type flux_fraction: float, optional
    :return: The model with the flux minimality constraint added.
    :rtype: cobra.Model
    """

    def recursive_sum(x):
        """Recursive version of sum.

        Sum operations are slow on variables. By doing it recursively, we can use
        2log(length(x)) additions instead of length(x).
        We could exceed the recursion depth, but it is not very likely with normal model sizes.
        """
        length = len(x)
        if length == 1:
            return x[0]
        elif length == 2:
            return x[0] + x[1]
        else:
            split = length // 2
            return recursive_sum(x[split:]) + recursive_sum(x[:split])

    # Calculate flux sum for the objective.
    flux_sum_value = cobra.flux_analysis.pfba(model, objective_fraction, objective)

    # Create a total flux variable (fw + rev fluxes, same as pFBA does it in Cobra.)
    variables = ((r.forward_variable, r.reverse_variable) for r in model.reactions)
    sum_variable = recursive_sum(list(itertools.chain.from_iterable(variables)))
    # Create and add the constraint to the model.
    constraint = model.problem.Constraint(
        expression=sum_variable,
        name="total_flux_sum_bound",
        lb=0,
        ub=flux_sum_value.objective_value * flux_fraction,
    )
    model.add_cons_vars([constraint])
    return model


def report_from_sample_frame(sample_frame, report_list, out_file=None):
    """take sampling output and summarize (mean) and report in dict form

    :param sample_frame:
    :param report_list:
    :param out_file:
    :return:
    """
    report_dict = {}
    if out_file:
        mean_frame = sample_frame.mean().to_frame()
        std_frame = sample_frame.std().to_frame()
        #combi_frame = mean_frame.merge(std_frame, left_index=True)
        combi_frame = pd.merge(mean_frame, std_frame, right_index=True, left_index=True)
        # print(mean_frame)
        # print(std_frame)
        # print(combi_frame)
        combi_frame.to_csv(out_file)
    for reaction in report_list:
        report_dict[reaction] = {}
        reaction_df = sample_frame[reaction]
        mean_value = reaction_df.mean()
        std_value = reaction_df.std()
        report_dict[reaction]['mean'] = mean_value
        report_dict[reaction]['sd'] = std_value
    return report_dict


def try_out_separate_pools(model_file):
    trial_link_reac_to_category_file = "../input_files/trial_proteome_file.tsv"
    cat_dict = {"ETC":0.01, "TCA":0.01, 'other':0.225}
    enz_model = get_enzyme_constraint_sep_protein_pools(trial_link_reac_to_category_file, model_file, cat_dict, verbose=True)
    enz_model.reactions.EX_pi_e.lower_bound = -0.14
    solution = enz_model.optimize()
    print(solution)

    #NB it will important to set kcat_MW values for the relevant reactions
    return None


def run_sampling(model_file, sampling_settings, out_file=None):
    enz_model = get_enzyme_constraint_model(model_file)
    enz_model.reactions.EX_pi_e.lower_bound = -0.14
    sample_frame = sample_model_with_pfba(enz_model, sampling_settings)
    result_dict = report_from_sample_frame(sample_frame, [], out_file=out_file)
    return result_dict


def json_write(path, dictionary):
    """Writes a JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path:   The path of the JSON file that shall be written
    * dictionary: The dictionary which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(dictionary, indent=4)
    with open(path, "w", encoding="utf-8") as f:
        f.write(json_output)


def export_ecmpy_model(model_object, output_json_file):
    """Take a cobra model of an ecmpy model and create an export in json format

    :param model_object:
    :param output_json_file:
    :return:
    """
    cobra.io.save_json_model(model_object, output_json_file)
    dictionary_model = json_load(output_json_file)
    for eachreaction in range(len(dictionary_model['reactions'])):
        reaction_id = dictionary_model['reactions'][eachreaction]['id']
        dictionary_model['reactions'][eachreaction]['kcat'] = enz_model.reactions.get_by_id(reaction_id).kcat
        dictionary_model['reactions'][eachreaction]['kcat_MW'] = enz_model.reactions.get_by_id(reaction_id).kcat_MW
    json_write(output_json_file, dictionary_model)
    return output_json_file


if __name__ == "__main__":
    strep_ec_model_file = "./model/iDG1237_updated_irr_enz_constraint_v4.json"
    non_constrained_model = "../Model_files/iDG1237_updated_with_uniprot.xml"

    strep_ec_model_file = "./model/iDG1237_updated_irr_enz_constraint_autopacman_v1.json"

    enz_model = get_enzyme_constraint_model(strep_ec_model_file)

    try_basic_enz_model_functions(strep_ec_model_file)

    #quick_try_non_constrained_model(non_constrained_model)

    #try_out_separate_pools(strep_ec_model_file)

    sampler_settings = {'objective': 'Growth',
                     'objective_fraction': 1,
                     'flux_fraction': 1.25,
                     'fraction_of_optimum': 0.99,
                     'n_samples': 1000}

    sample_file = 'sampling_pFBA_test_autopacman_v1_phos_0_14_frac_opt_0_99.csv'
    #run_sampling(strep_ec_model_file, sampler_settings, sample_file)
