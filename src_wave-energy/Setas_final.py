import matplotlib.gridspec as gridspec
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Carregar as datas de eventos de ressaca de um arquivo CSV
datas_ressaca = pd.read_csv('/p1-nemo/rtecchio/teste_chico/dias_ressaca.csv', sep=';')
# Combinar as colunas "ano", "mes" e "dia" para formar a data completa
datas_ressaca["data"] = pd.to_datetime(datas_ressaca["ano"].astype(str) + "-" + datas_ressaca["mes"].astype(str) + "-" + datas_ressaca["dia"].astype(str))
# Selecionar apenas a coluna com as datas completas
datas_ressaca = datas_ressaca["data"]

# Carregar os dados diários de um arquivo CSV
dados = pd.read_csv('/p1-nemo/rtecchio/teste_chico/variaveis_branco.csv', sep=',', parse_dates=["time"])

# Selecionar os dias concomitantes, dois dias antes e um dia depois das datas de eventos de ressaca
datas_selecionadas = []
for data_ressaca in datas_ressaca:
    dias_selecionados = dados[(dados["time"] >= data_ressaca - pd.Timedelta(days=2)) &
                              (dados["time"] <= data_ressaca + pd.Timedelta(days=1))]
    datas_selecionadas.extend(dias_selecionados["time"])

# Remover duplicatas das datas selecionadas
datas_selecionadas = list(set(datas_selecionadas))

# Ordenar as datas selecionadas
datas_selecionadas.sort()

# Plotando a série temporal de 2010 a 2019
# Definir as cores e rótulos para cada ponto
cores = ['r', 'g', 'b']
rotulos = ['Ponto 1', 'Ponto 2', 'Ponto 3']


# Percorrer cada ponto e criar uma imagem separada
for i, ponto in enumerate([1, 2, 3]):
    # Filtrar os dados para o ponto atual
    ponto_dados = dados[dados["ponto"] == ponto]
    
    # Configurar o gráfico
    plt.figure()
    plt.xlabel("Data")
    plt.ylabel("P")
    plt.title(f"Série Temporal de P - Ponto {ponto} (2010-2019)")
    
    # Plotar a série temporal de PPer
    plt.plot(ponto_dados["time"], ponto_dados["P"], color=cores[i], label=rotulos[i])
    
    # Adicionar as barras coloridas para identificar as datas dos eventos de ressaca
    for data_ressaca in datas_ressaca:
        plt.axvspan(data_ressaca - pd.Timedelta(days=1), data_ressaca + pd.Timedelta(days=3), facecolor='black', alpha=0.5)
    
    # Adicionar a legenda
   # plt.legend()
    
    # Salvar a imagem em um arquivo separado
    plt.savefig(f"SerieTemporal_P_Ponto{ponto}.png")
    
    # Exibir o gráfico
    plt.show()



# Loop através de cada data de evento de ressaca
for data_ressaca in datas_ressaca:
    # Filtrar os dados para o evento de ressaca atual
    evento_ressaca = dados[(dados["time"] >= data_ressaca - pd.Timedelta(days=2)) &
                           (dados["time"] <= data_ressaca + pd.Timedelta(days=1))]
    
    # Plotar as curvas para cada ponto com cor e rótulo distintos
    for i, ponto in enumerate([1, 2, 3]):
        ponto_ressaca = evento_ressaca[evento_ressaca["ponto"] == ponto]
        plt.plot(ponto_ressaca["time"], ponto_ressaca["PPer"], color=cores[i], label=rotulos[i])
    
    plt.xlabel("Data")
    plt.ylabel("PPer (kW/m)")
    plt.title(f"Evento de Ressaca - {data_ressaca.date()}")
    plt.legend()
    
    # Salvar o gráfico em um arquivo
    nome_arquivo = f"evento_ressaca_PPer{str(data_ressaca.date())}.png"
    plt.savefig(nome_arquivo)
    plt.show()
    
    # Limpar o gráfico para o próximo evento
    plt.clf()


# Plotando setas para a energia nos pontos

import xarray as xr
import pandas as pd
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature, COASTLINE, BORDERS
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cmocean as cmo

# Profundidades para deixar a figura mais bonita
ds = xr.open_dataset('/p1-nemo/rtecchio/Dados/GEBCO/gebco_2023_costa_s_se.nc')
prof2 = ds['elevation'][:]
prof_filtred = prof2.where(prof2 <= 0)
prof_filtred_min = prof_filtred.min().values
prof_filtred_max = prof_filtred.max().values

# Pontos
tan_ptos = pd.read_csv('/p1-nemo/rtecchio/teste_chico/pontos_branco_dir_perp.csv', sep=';', decimal=',')
lat_pto = tan_ptos['lat']
lon_pto = tan_ptos['lon']
dir_perp = tan_ptos['DIR_NORM']
dir_par= tan_ptos['DIR_PAR']

# Converter a direção para radianos
#decompondo a direção
direcao_graus = dir_perp
direcao_graus2 = dir_par
mag = 2.5
# Converter a direção para radianos
direcao_radianos = np.deg2rad(direcao_graus)
direcao_radianos2 = np.deg2rad(direcao_graus2)
# Calcular as componentes x e y do vetor
u = mag * np.cos(direcao_radianos)
v = mag * np.sin(direcao_radianos)
#PPar
u2 = mag * np.cos(direcao_radianos2)
v2 = mag * np.sin(direcao_radianos2)

# Definir o número de dias antes e depois do evento de ressaca
dias_antes = 2
dias_depois = 1


def map_features(ax):
    ax.add_feature(COASTLINE)
    ax.add_feature(BORDERS, edgecolor='#383838')
    return ax

def Brazil_states(ax):    
    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                  name='admin_1_states_provinces_lines')
    _ = ax.add_feature(states, edgecolor='#383838')
    
    cities = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                  name='populated_places')
    _ = ax.add_feature(cities)
    
def grid_labels_params(ax,i):
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5,linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    if i not in [0,3]:
        gl.left_labels = False
    gl.xlabel_style = {'size': 12, 'color': '#383838'}
    gl.ylabel_style = {'size': 12, 'color': '#383838'}
    ax.spines['geo'].set_edgecolor('#383838')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return ax

# Iterar sobre as datas de eventos de ressaca
for data_ressaca in datas_ressaca:
    # Calcular a janela de tempo para o evento de ressaca
    inicio_janela = data_ressaca - pd.Timedelta(days=dias_antes)
    fim_janela = data_ressaca + pd.Timedelta(days=dias_depois)
    
    # Filtrar os dados dentro da janela de tempo
    dados_janela = dados[(dados["time"] >= inicio_janela) & (dados["time"] <= fim_janela)]
    
    # Calcular a média diária dos dados para cada ponto
    dados_diarios = dados_janela.groupby([dados_janela["time"].dt.date, dados_janela["ponto"]]).mean()
    
    # Obter os valores de PPer e PPar para cada ponto
    PPer_ponto1 = dados_diarios.loc[(slice(None), 1), "PPer"].values
    PPer_ponto2 = dados_diarios.loc[(slice(None), 2), "PPer"].values
    PPer_ponto3 = dados_diarios.loc[(slice(None), 3), "PPer"].values
    
    PPar_ponto1 = dados_diarios.loc[(slice(None), 1), "Ppar"].values
    PPar_ponto2 = dados_diarios.loc[(slice(None), 2), "Ppar"].values
    PPar_ponto3 = dados_diarios.loc[(slice(None), 3), "Ppar"].values
    
    # Criar uma figura para o evento de ressaca
    fig = plt.figure(figsize=(9, 8))
    
    # Criar uma grade de subplots
    gs = gridspec.GridSpec(2, 2, figure=fig, wspace=-0.2, hspace=0.2)
    
    # Criar subplots para os quatro dias
    for i, dia in enumerate(range(-dias_antes, dias_depois + 1)):
        #Incrementar 1 ao valor de dia para o título
        dia += 2
        ax = fig.add_subplot(gs[i], projection=ccrs.PlateCarree())
        
        # Configurar o mapa e suas características
        ax.set_extent([-42, -37, -22, -17.7])
        map_features(ax)
        Brazil_states(ax)
        grid_labels_params(ax, i)
        cf1 = ax.contourf(prof_filtred.lon[::20], prof_filtred.lat[::20], prof_filtred[::20, ::20], cmap='cmo.deep_r')
        
        ax.contourf(prof_filtred.lon[::20], prof_filtred.lat[::20], prof_filtred[::20, ::20], cmap='cmo.deep_r')
        
        # Plotar as setas de energia
        lon_pto1 = list(map(float, lon_pto.values))
        lat_pto1 = list(map(float, lat_pto.values))
        
        color_energy = PPer_ponto1 if i in [0, 1] else PPer_ponto2 if i == 2 else PPer_ponto3
        color_energy2 = PPar_ponto1 if i in [0, 1] else PPar_ponto2 if i == 2 else PPar_ponto3
        
        cmap_dia = cm.get_cmap('jet')
        norm_dia = colors.TwoSlopeNorm(vmin=-10, vcenter=0, vmax=50)
        colors_dia = cmap_dia(norm_dia(color_energy))
        colors_dia2 = cmap_dia(norm_dia(color_energy2))
        
        ax.quiver(lon_pto1, lat_pto1,  u, v, color=colors_dia, edgecolor='k', linewidth=1, pivot='tip', scale=20, width=0.015)
        ax.quiver(lon_pto1, lat_pto1,  u2, v2, color=colors_dia2, edgecolor='k', linewidth=1, pivot='tip', scale=20, width=0.015)
        
        # Adicionar título ao subplot com a data do dia
        titulo = f'Dia: {inicio_janela + pd.Timedelta(days=dia):%Y-%m-%d}'
        ax.set_title(titulo, fontsize=12)
    
    # Adicionar barra de cores para PPer
    cbar_PPer = fig.add_axes([0.057, 0.05, 0.40, 0.02])  # Posição da colorbar PPer
    cbper = fig.colorbar(cm.ScalarMappable(norm=norm_dia, cmap=cmap_dia), cax=cbar_PPer, orientation='horizontal')
    cbper.ax.tick_params(labelsize=8)
    cbper.set_label('PPer', labelpad=0, fontsize=12, fontweight='bold')
    
    # Adicionar barra de cores para PPar
    cbar_PPar = fig.add_axes([0.57, 0.05, 0.40, 0.02])  # Posição da colorbar PPar
    cbpar = fig.colorbar(cm.ScalarMappable(norm=norm_dia, cmap=cmap_dia), cax=cbar_PPar, orientation='horizontal')
    cbpar.ax.tick_params(labelsize=8)
    cbpar.set_label('PPar', labelpad=0, fontsize=12, fontweight='bold')
    
    plt.subplots_adjust(left=0.02, right=1, bottom=0.15, top=0.90)
    
    # Salvar a figura em um arquivo separado com o nome baseado na data do evento de ressaca
    nome_arquivo = f'P_setas_grid2_{data_ressaca.strftime("%Y%m%d")}.png'
    plt.savefig(nome_arquivo, dpi=300)
    
    # Fechar a figura para liberar memória
    plt.close(fig)
