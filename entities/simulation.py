from entities.chemicalProcess import ChemicalProcess

class Simulation:

    def __init__(self,input,sys_configs,process_configs):
        self.problem_inputs = self.validate_inputs(input)
        ##Sys Configs
        self.rec_stream_initial_guess = sys_configs['rec_stream_initial_guess']
        self.rec_compositions_initial_guess = [sys_configs['rec_compositions_initial_guess']]*process_configs['N_components']
        # self.rec_compositions_initial_guess = sys_configs['rec_compositions_initial_guess']
        self.max_iterations = sys_configs['max_iterations']
        self.convergence_threshold = sys_configs['convergence_threshold']
        ##Inputs
        self.Fo = input['Fo']
        self.Win = [input['Xoa'],input['Xob'],input['Xoc'],input['Xod']]
        self.Pr = self.bar_to_pascal(input['Pr'])
        self.Tr = input['Tr']
        self.Pf = self.bar_to_pascal(input['Pf'])
        self.Tf = input['Tf']
        self.Cs = input['Cs']
        ##Process Configs
        self.Vr = process_configs['Vr']
        self.Kor = process_configs['Kor']
        self.Ea = process_configs['Ea']
        self.reaction_coefficients = process_configs['reaction_coefficients']
        self.reaction_model = process_configs['reaction_model']
        self.lv_equilibrum_model = process_configs['lv_equilibrum_model']
        self.antoinecoefficients = process_configs['antoinecoefficients']
        ##Composed -- to be changed
        self.lv_methodinput = [self.Tf,self.antoinecoefficients]

    def bar_to_pascal(self,P):
        return P*(10**5)

    def validate_inputs(self,input):
        problem_inputs = []
        if input['Fo'] > 200 or input['Fo'] < 50: problem_inputs.append("Fo inserted out of avaliable range")
        if input['Pr'] > 14 or input['Pr'] < 10: problem_inputs.append("Pr inserted out of avaliable range")
        if input['Tr'] > 1250 or input['Tr'] < 850: problem_inputs.append("Tr inserted out of avaliable range")
        if input['Pf'] > 7 or input['Pf'] < 3: problem_inputs.append("Pf inserted out of avaliable range")
        if input['Tf'] > 700 or input['Tf'] < 300: problem_inputs.append("Tf inserted out of avaliable range")
        if input['Cs'] > 0.8 or input['Cs'] < 0: problem_inputs.append("Cs inserted out of avaliable range")
        if input['Xoa'] > 1 or input['Xoa'] < 0: problem_inputs.append("Xoa inserted out of avaliable range")
        if input['Xob'] > 1 or input['Xob'] < 0: problem_inputs.append("Xob inserted out of avaliable range")
        if input['Xoc'] > 1 or input['Xoc'] < 0: problem_inputs.append("Xoc inserted out of avaliable range")
        if input['Xod'] > 1 or input['Xod'] < 0: problem_inputs.append("Xod inserted out of avaliable range")
        if round(input['Xoa']+input['Xob']+input['Xoc']+input['Xod'],4) != 1.0000: problem_inputs.append("Molar ratios do not sum zero. Check compositions inserted.")
        return problem_inputs

    def calculate_results(self):
        N_iteration=0
        while self.max_iterations > N_iteration:
            simul = ChemicalProcess(self.rec_stream_initial_guess,self.rec_compositions_initial_guess)
            simul.evaluate(self.Fo,self.Win,self.Vr,self.Pr,self.Tr,
                            self.reaction_coefficients,self.reaction_model,self.Kor,self.Ea,
                            self.lv_equilibrum_model,self.Pf,self.lv_methodinput,self.Cs)
            self.rec_stream_initial_guess = simul.F[6]
            self.rec_compositions_initial_guess = simul.W[6]
            if simul.residual < self.convergence_threshold:
                print('Rodou legal')
                return simul
            N_iteration = N_iteration+1
        return None
    
    def run_simulation(self):
        if len(self.problem_inputs) > 0:
            self.write_warning()
            return None
        last_iteration = self.calculate_results()
        self.write_output(last_iteration)

    def write_warning(self):
        f = open("warning.txt", "w")
        for warning in self.problem_inputs:
            f.write(warning + "\n")
        f.close()

    def write_output(self,last_iteration_data):
        f = open("output.txt", "w")
        if last_iteration_data == None:
            f.write("Calculation did not converge.")
        else:
            data_to_write = self.format_output(last_iteration_data)
            f.write(data_to_write)
        f.close()

    def format_output(self,last_iteration_data):
        template = self.get_output_template()
        output_text = self.fill_output_text(template,last_iteration_data)
        return output_text

    def get_output_template(self):
        f2 = open("./configs/output_config.txt", "r")
        output_template = f2.read()
        return output_template

    def fill_output_text(self,output_text,last_iteration_data):
        for i in range(len(last_iteration_data.F)):
            output_text = output_text.replace(f"F{i}",self.format_result_numbers(last_iteration_data.F[i]))
            for j in range(len(last_iteration_data.W[0])):
                output_text = output_text.replace(f"X{i}{j}",self.format_result_numbers(last_iteration_data.W[i][j]))
        return output_text

    def format_result_numbers(self,x):
        return str(round(x,3))
