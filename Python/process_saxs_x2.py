#!/usr/bin/env python
# coding: utf-8

# In[100]:


###########################################################################
#           Code to process SAXS output files to CSV files                #
#                       and calculate the X^2                             #
#                    Author: Renato D. Cunha, PhD                         #
#                     E-mail: renatodias@ub.edu                           #
#                          Version: 1.0                                   #
###########################################################################
#                                                                         #
# READ ME:                                                                #
# 1 - All the .abs files need to be in the same folder                    #
# 2 - The data is organized by default with ',' as separator              #
# 3 - Feel free to adapt the code to your needs                           #
#                                                                         #
###########################################################################


# In[1]:


import os
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize


# In[2]:


def converter_arquivos_para_csv():
    # Obtém o diretório atual
    pasta_entrada = os.getcwd()

    # Verifica se o diretório de saída existe, senão cria
    pasta_saida = 'csv_output'
    os.makedirs(pasta_saida, exist_ok=True)

    # Lista todos os arquivos com extensão .abs no diretório de entrada
    arquivos_abs = [f for f in os.listdir(pasta_entrada) if f.endswith('.abs')]

    for arquivo_abs in arquivos_abs:
        # Constrói o caminho completo para o arquivo de entrada
        caminho_abs = os.path.join(pasta_entrada, arquivo_abs)

        # Lista para armazenar linhas processadas
        linhas_processadas = []

        # Flag para indicar se é a primeira linha
        primeira_linha = True

        # Lê o arquivo .abs linha a linha
        with open(caminho_abs, 'r') as file:
            for linha in file:
                # Pula a primeira linha que contém os títulos das colunas
                if primeira_linha:
                    primeira_linha = False
                    continue

                # Extrai os valores manualmente com base na posição
                valor1 = float(linha[:21].replace('D', 'E').strip())  # Substitui 'D' por 'E' na notação científica
                valor2 = float(linha[21:42].replace('D', 'E').strip())  # Substitui 'D' por 'E' na notação científica

                # Adiciona os valores processados à lista
                linhas_processadas.append([valor1, valor2])

        # Converte a lista de linhas processadas para DataFrame
        dados = pd.DataFrame(linhas_processadas, columns=['Coluna1', 'Coluna2'])

        # Constrói o caminho completo para o arquivo de saída (CSV)
        caminho_csv = os.path.join(pasta_saida, f"{os.path.splitext(arquivo_abs)[0]}.csv")

        # Salva os dados no formato CSV com vírgula como separador
        dados.to_csv(caminho_csv, index=False, sep=',')

        print(f"Conversão concluída para {caminho_csv}")


# In[3]:


def calcular_media_colunas_csv():
    # Obtém o diretório atual
    pasta_trabalho = os.path.join(os.getcwd(), 'csv_output')

    # Lista todos os arquivos com extensão .csv na pasta "csv_output"
    arquivos_csv = [f for f in os.listdir(pasta_trabalho) if f.endswith('.csv')]

    # Lista para armazenar os DataFrames de cada arquivo
    dfs = []

    # Loop sobre os arquivos CSV
    for arquivo_csv in arquivos_csv:
        # Constrói o caminho completo para o arquivo CSV
        caminho_csv = os.path.join(pasta_trabalho, arquivo_csv)

        # Lê os dados do arquivo CSV
        dados = pd.read_csv(caminho_csv, index_col=0)  # Assumindo que a primeira coluna é o índice

        # Adiciona o DataFrame à lista
        dfs.append(dados)

        # Print para mostrar o arquivo que está sendo usado
        print(f"Usando arquivo: {arquivo_csv}")

    # Calcula a média final para todas as colunas através dos arquivos
    media_final_colunas = pd.concat(dfs, axis=1).mean(axis=1)

    # Cria um DataFrame com a média final das colunas
    df_resultado = pd.DataFrame({'Média Colunas': media_final_colunas})

    # Print para visualizar o DataFrame
    print("\nDataFrame Resultado:")
    print(df_resultado)

    # Salva o DataFrame no arquivo Resultado_Medias.csv
    caminho_saida = os.path.join(pasta_trabalho, 'SAXS_avg.csv')
    df_resultado.to_csv(caminho_saida)


# In[4]:


def processar_e_escalar_saxs_avg():
    # Obtém o diretório atual
    pasta_trabalho = os.getcwd()

    # Caminho completo para o arquivo SAXS_avg.csv na pasta "csv_output"
    caminho_saxs_avg = os.path.join(pasta_trabalho, 'csv_output', 'SAXS_avg.csv')

    # Verifica se o arquivo SAXS_avg.csv existe
    if os.path.exists(caminho_saxs_avg):
        # Carrega o DataFrame do arquivo SAXS_avg.csv
        df_saxs_avg = pd.read_csv(caminho_saxs_avg)

        # Escala os valores da coluna 1 por 10 e da coluna 2 divide por 3
        df_saxs_avg.iloc[:, 0] *= 10
        df_saxs_avg.iloc[:, 1] /= 2.206

        # Salva o DataFrame no novo arquivo SAXS_avg_processado.csv
        caminho_novo_saxs_avg = os.path.join(pasta_trabalho, 'SAXS_avg_processado.csv')
        df_saxs_avg.to_csv(caminho_novo_saxs_avg, index=False, header=['q_theo', 'I_theo'])

        # Print para verificar o DataFrame processado
        print("\nDataFrame Processado:")
        print(df_saxs_avg)

        print(f"\nNovo arquivo SAXS_avg_processado.csv salvo em: {caminho_novo_saxs_avg}")
    else:
        print("O arquivo SAXS_avg.csv não foi encontrado na pasta 'csv_output'.")



# In[5]:


def interpolar_e_salvar_saxs_theo(arquivo_exp, arquivo_avg_processado):
    # Carregar os dois arquivos, ignorando a primeira linha
    df_exp = pd.read_csv(arquivo_exp, skiprows=1)
    df_avg_processado = pd.read_csv(arquivo_avg_processado, skiprows=1)

    # Extrair valores de x e y dos DataFrames
    x_exp, y_exp = df_exp.iloc[:, 0], df_exp.iloc[:, 1]
    x_avg_processado, y_avg_processado = df_avg_processado.iloc[:, 0], df_avg_processado.iloc[:, 1]

    # Interpolar os valores de y_avg_processado para os pontos de x_exp
    interpolador = interp1d(x_avg_processado, y_avg_processado, kind='linear', fill_value='extrapolate')
    y_theo_interpolado = interpolador(x_exp)

    # Criar DataFrame com os valores interpolados
    df_theo_interpolado = pd.DataFrame({'q_theo': x_exp, 'I_theo_interpolado': y_theo_interpolado})

    # Salvar o DataFrame no arquivo SAXS_theo.csv
    caminho_saxs_theo = os.path.join(os.getcwd(), 'SAXS_theo.csv')
    df_theo_interpolado.to_csv(caminho_saxs_theo, index=False)

    print(f"\nArquivo SAXS_theo.csv salvo com os valores interpolados.")


# In[6]:


def calcular_x2_e_salvar(arquivo_exp, arquivo_theo):
    # Carregar os dois arquivos, ignorando a primeira linha
    df_exp = pd.read_csv(arquivo_exp, skiprows=1)
    df_theo = pd.read_csv(arquivo_theo, skiprows=1)

    # Garantir que ambos os DataFrames tenham o mesmo número de linhas
    min_linhas = min(df_exp.shape[0], df_theo.shape[0])
    df_exp = df_exp.head(min_linhas)
    df_theo = df_theo.head(min_linhas)

    # Garantir que SAXS_theo.csv tenha 2 colunas e SAXS_exp_fixed.csv tenha 3 colunas
    df_theo = df_theo.iloc[:, :2]
    df_exp = df_exp.iloc[:, :3]

    # Double-check: garantir que a primeira coluna dos dois arquivos seja igual
    if not df_exp.iloc[:, 0].equals(df_theo.iloc[:, 0]):
        raise ValueError("As primeiras colunas dos arquivos não são iguais.")

    # Atribuir variáveis para facilitar o cálculo
    I_theo = df_theo.iloc[:, 1]
    I_exp = df_exp.iloc[:, 1]
    I_error = df_exp.iloc[:, 2]

    # Calcular x_2 para todos os pontos
    x_2 = (((I_exp - (0.682103 * I_theo + 0.000009) / I_error) ** 2

    # Criar DataFrame com os valores calculados
    df_x2_calc = pd.DataFrame({'q_exp': df_exp.iloc[:, 0], 'x_2': x_2})

    # Salvar o DataFrame no arquivo x2_calc.csv
    caminho_x2_calc = os.path.join(os.getcwd(), 'x2_calc.csv')
    df_x2_calc.to_csv(caminho_x2_calc, index=False)

    # Printar os valores de x_2
    print("\nValores de x_2 para cada ponto:")
    print(df_x2_calc)

    # Calcular e printar a média de todos os x_2 calculados
    media_x2 = x_2.sum() / 180
    print(f"\nMédia de x_2: {media_x2}")

    print(f"\nArquivo x2_calc.csv salvo com os valores de x_2.")


# In[16]:


def calcular_x2_e_salvar_otimizado(arquivo_exp, arquivo_theo):
    # Carregar os dois arquivos, ignorando a primeira linha
    df_exp = pd.read_csv(arquivo_exp, skiprows=1)
    df_theo = pd.read_csv(arquivo_theo, skiprows=1)

    # Garantir que ambos os DataFrames tenham o mesmo número de linhas
    min_linhas = min(df_exp.shape[0], df_theo.shape[0])
    df_exp = df_exp.head(min_linhas)
    df_theo = df_theo.head(min_linhas)

    # Garantir que SAXS_theo.csv tenha 2 colunas e SAXS_exp_fixed.csv tenha 3 colunas
    df_theo = df_theo.iloc[:, :2]
    df_exp = df_exp.iloc[:, :3]

    # Double-check: garantir que a primeira coluna dos dois arquivos seja igual
    if not df_exp.iloc[:, 0].equals(df_theo.iloc[:, 0]):
        raise ValueError("As primeiras colunas dos arquivos não são iguais.")

    # Atribuir variáveis para facilitar o cálculo
    q_exp = df_exp.iloc[:, 0]
    I_theo = df_theo.iloc[:, 1]
    I_exp = df_exp.iloc[:, 1]
    I_error = df_exp.iloc[:, 2]

    # Função objetivo para minimizar
    def objetivo(params):
        c, b = params
        x_2 = ((I_exp - (c * I_theo + b))/ I_error) ** 2
        return x_2.sum()
                
    # Valor inicial para c e b
    parametros_iniciais = [1.0, 1.0]

    # Otimização
    resultado_otimizacao = minimize(objetivo, parametros_iniciais, method='BFGS')

    # Parâmetros otimizados
    c_otimizado, b_otimizado = resultado_otimizacao.x

    print(f"\nParâmetros otimizados: c = {c_otimizado:.4f}, b = {b_otimizado:.4f}")

    # Calcular x_2 para todos os pontos usando os parâmetros otimizados
    x_2_otimizado = (((I_exp - ((c_otimizado * I_theo) + b_otimizado))) / I_error) ** 2

    # Criar DataFrame com os valores calculados e os parâmetros otimizados
    df_x2_calc = pd.DataFrame({'q_exp': q_exp, 'x_2': x_2_otimizado, 'c_otimizado': c_otimizado, 'b_otimizado': b_otimizado})

    # Salvar o DataFrame no arquivo x2_calc.csv
    caminho_x2_calc = os.path.join(os.getcwd(), 'x2_calc.csv')
    df_x2_calc.to_csv(caminho_x2_calc, index=False)

    # Printar os valores de x_2 otimizado e os parâmetros otimizados
    print("\nValores de x_2 otimizado para cada ponto e parâmetros otimizados:")
    print(df_x2_calc)

    # Calcular e printar a média de todos os x_2 otimizados
    media_x2_otimizado = x_2_otimizado.mean()
    
    print(f"\nMédia de x_2 otimizado: {media_x2_otimizado}")

    print(f"\nArquivo x2_calc.csv salvo com os valores de x_2 otimizado e parâmetros otimizados.")


# In[10]:





# In[ ]:





# In[17]:


# Exemplo de uso
converter_arquivos_para_csv()
calcular_media_colunas_csv()
processar_e_escalar_saxs_avg()
interpolar_e_salvar_saxs_theo('SAXS_exp.csv', 'SAXS_avg_processado.csv')
calcular_x2_e_salvar_otimizado('SAXS_exp_fixed.csv', 'SAXS_theo.csv')
#calcular_x2_e_salvar('SAXS_exp_fixed.csv', 'SAXS_theo.csv')


# In[ ]:





# In[ ]:




