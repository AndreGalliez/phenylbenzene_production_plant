from entities.chemicalProcess import ChemicalProcess

class Simulation:
    """
        Responsável por determinar a ordem em que os equipamentos são calculados, 
        chamar o calculo e passar adiante os outputs de equipamentos que são inputs de outros.
        Argumentos:
            input (dict): Input com os dados de entrada do processo fornecidos pelo usuário. 
            sys_configs (dict): Dicionário contendo configurações do processo de cálculo:
                (Quantas iterações permitir, qual chute inicial usar para a corrente de reciclo, qual diferença tolerar para convergência)
            process_configs (dict): Dicionário contendo configurações do processo químico:
                (Constante padrão reacional Ko, Energia de ativação Ea, Coeficientes de Reação e Parâmetros dos modelos de equilibrio LV).
        Atributos:
        rec_stream_initial_guess = Chute inicial da vazão de reciclo (configuração de cálculo).
        rec_compositions_initial_guess = Chute inicial das composições de reciclo (configuração de cálculo).
        self.max_iterations = Número máximo de iterações permitidas (configuração de cálculo).
        self.convergence_threshold = Critério limite de convergência (configuração de cálculo).
        Fo = (float) Vazão de entrada (input).
        Win = (list(float)) Composições de entrada (input).
        Pr = (float) Pressão no reator (input).
        Tr = (float) Temperatura no reator (input).
        Pf = (float) Pressão no flash (input).
        Tf = (float) Temperatura no flash (input).
        Cs = (float) Razão de reciclo (input).
        Vr = (float) Volume do reator (configuração de processo).
        Kor = list(list((float))) Constantes padrão de reação (configuração de processo).
        Ea = list(list((float))) Energia de ativação (configuração de processo).
        reaction_coefficients list(list((float))) Coeficientes reacionais (configuração de processo).
        elv_coefficients = list(list((float))) Parâmetros dos modelos de equilibrio LV (configuração de processo).
        Métodos:
            bar_to_pascal()
                : Converte pressões em bar (dado de entrada) para Pascal (usado no calculo).
            validate_inputs()
                : Realiza a validação dos inputs, garantindo que eles estajam dentro das faixas permitdas.
            calculate_results()
                : Instancia os objetos ChemicalProcess e realiza os cálculos de processo de forma iterativa até a convergência.
            write_warning()
                : Escreve um arquivo com avisos de inputs incorretos. Utilizado com os dados de entrada do usuário não são apropriados.
            format_result_numbers()
                : Formata os resultados obtidos para 3 casas decimais.
            get_output_template()
                : Lê o arquivo que define o formato a ser usado no arquivo de saída.
            fill_output_text()
                : Preenche o texto do arquivo de saída com as informações obtidas dos cálculos.
            format_output()
                : Wrapper que lê o template com o formato do arquivo de saída, preenche com os dados de processo e retorna o conteúdo completo. 
            write_output()
                : Escreve o conteúdo do arquivo de saída em um arquivo de texto.
            run_simulation()
                : Wrapper que chama a simulação, roda todas as outras funções na ordem correta e escreve os arquivos pertinentes.
    """

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
        self.elv_coefficients = process_configs['elv_coefficients']

    def bar_to_pascal(self,P):
        return P*(10**5)

    def validate_inputs(self,input):
        """
            Realiza a validação dos inputs, garantindo que eles estajam dentro das faixas permitdas.
            Argumentos:
                input (dict): Input com os dados de entrada do processo fornecidos pelo usuário.
            Retorna:
                problem_inputs (list): list com os avisos a serem escritos em um arquivo "warning" caso haja algum.
        """
        problem_inputs = []
        if input['Fo'] > 200 or input['Fo'] < 50: problem_inputs.append("Fo inserted out of allowed range.")
        if input['Pr'] > 14 or input['Pr'] < 10: problem_inputs.append("Pr inserted out of allowed range.")
        if input['Tr'] > 1250 or input['Tr'] < 850: problem_inputs.append("Tr inserted out of allowed range.")
        if input['Pf'] > 12 or input['Pf'] < 3: problem_inputs.append("Pf inserted out of allowed range.")
        if input['Tf'] > 700 or input['Tf'] < 300: problem_inputs.append("Tf inserted out of allowed range.")
        if input['Cs'] > 0.8 or input['Cs'] < 0: problem_inputs.append("Cs inserted out of allowed range.")
        if input['Xoa'] > 1 or input['Xoa'] < 0: problem_inputs.append("Xoa inserted out of allowed range.")
        if input['Xob'] > 1 or input['Xob'] < 0: problem_inputs.append("Xob inserted out of allowed range.")
        if input['Xoc'] > 1 or input['Xoc'] < 0: problem_inputs.append("Xoc inserted out of allowed range.")
        if input['Xod'] > 1 or input['Xod'] < 0: problem_inputs.append("Xod inserted out of allowed range.")
        if round(input['Xoa']+input['Xob']+input['Xoc']+input['Xod'],4) != 1.0000: problem_inputs.append("Molar ratios do not sum zero. Check compositions inserted.")
        return problem_inputs

    def calculate_results(self):
        """
            Instancia os objetos ChemicalProcess e realiza os cálculos de processo de forma iterativa até a convergência.
            Argumentos:
            Retorna:
                Objeto com a iteração do processo que convergiu, ou None caso não ocorra convergência (ChemicalProcess)
        """
        N_iteration=0
        while self.max_iterations > N_iteration:
            simul = ChemicalProcess(self.rec_stream_initial_guess,self.rec_compositions_initial_guess)
            simul.evaluate(self.Fo,self.Win,self.Vr,self.Pr,self.Tr,self.reaction_coefficients,self.Kor,self.Ea,
                            self.Pf,self.Tf,self.elv_coefficients,self.Cs)
            self.rec_stream_initial_guess = simul.F[6]
            self.rec_compositions_initial_guess = simul.W[6]
            if simul.residual < self.convergence_threshold:
                return simul
            N_iteration = N_iteration+1
        return None

    def write_warning(self):
        f = open("warning.txt", "w")
        for warning in self.problem_inputs:
            f.write(warning + "\n")
        f.close()

    def format_result_numbers(self,x):
        """
            Formata os resultados obtidos para 3 casas decimais.
            Argumentos:
                x (float): Numero a ser formatado
            Retorna:
                (float) numero formatado com 3 casas decimais
        """
        return str(round(x,3))

    def get_output_template(self):
        f2 = open("./configs/output_config.txt", "r")
        output_template = f2.read()
        return output_template

    def fill_output_text(self,output_text,last_iteration_data):
        """
            Preenche o texto do arquivo de saída com as informações obtidas dos cálculos.
            Argumentos:
                output_text (str): texto com o template do output a ser preenchido
                last_iteration_data (ChemicalProcess): objeto com a iteração que convergiu contendo as informações calculadas do processo.
            Retorna:
                (str) texto do arquivo de saída preenchido com as informações do cálculo
        """
        for i in range(len(last_iteration_data.F)):
            output_text = output_text.replace(f"F{i}",self.format_result_numbers(last_iteration_data.F[i]))
            for j in range(len(last_iteration_data.W[0])):
                output_text = output_text.replace(f"X{i}{j}",self.format_result_numbers(last_iteration_data.W[i][j]))
        return output_text

    def format_output(self,last_iteration_data):
        """
            Wrapper que lê o template com o formato do arquivo de saída, preenche com os dados de processo e retorna o conteúdo completo. 
            Argumentos:
                last_iteration_data (ChemicalProcess): objeto com a iteração que convergiu contendo as informações calculadas do processo.
            Retorna:
                (str) texto do arquivo de saída preenchido com as informações do cálculo
        """
        template = self.get_output_template()
        output_text = self.fill_output_text(template,last_iteration_data)
        return output_text

    def write_output(self,last_iteration_data):
        """
            Escreve o conteúdo do arquivo de saída em um arquivo de texto. 
            Implementa a lógica que lida a escrita de arquivos onde o cálculo não convergiu.
            Argumentos:
                last_iteration_data (ChemicalProcess): objeto com a iteração que convergiu contendo as informações calculadas do processo.
        """
        f = open("output.txt", "w")
        if last_iteration_data == None:
            f.write("Calculation did not converge.")
        else:
            data_to_write = self.format_output(last_iteration_data)
            f.write(data_to_write)
        f.close()

    def run_simulation(self):
        if len(self.problem_inputs) > 0:
            self.write_warning()
            return None
        last_iteration = self.calculate_results()
        self.write_output(last_iteration)

