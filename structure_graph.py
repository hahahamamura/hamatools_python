#!/usr/bin/env python3
"""
Script para gerar gráfico de barras a partir do output do STRUCTURE
Similar ao formato do CLUMPAK
Com cores consistentes baseadas na ancestralidade predominante das primeiras populações
"""

import re
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from matplotlib.patches import Rectangle
import argparse

def parse_structure_output(filename):
    """
    Parse do arquivo de output do STRUCTURE
    Retorna os dados de ancestralidade dos indivíduos e a ordem das populações
    """
    individuals = []
    populations = []
    clusters = []
    reading_individuals = False
    reading_given_pops = False
    num_clusters = None
    given_pop_order = []
    pop_info = {}
    
    with open(filename, 'r') as f:
        for line in f:
            # Detectar número de clusters
            if 'populations assumed' in line:
                num_clusters = int(re.search(r'(\d+)\s+populations assumed', line).group(1))
            
            # Detectar início da tabela de Given Pops
            if 'Given    Inferred Clusters' in line:
                reading_given_pops = True
                next(f)  # Pular linha "Pop    1    2    3..."
                continue
            
            # Ler tabela de Given Pops
            if reading_given_pops:
                line = line.strip()
                if line.startswith('--'):
                    reading_given_pops = False
                    continue
                
                # Parse: "  1:     0.004  0.984  0.008  0.004      504"
                match = re.match(r'\s*(\d+):\s+([\d.]+.*?)\s+(\d+)\s*$', line)
                if match:
                    pop_num = match.group(1)
                    n_individuals = int(match.group(3))
                    given_pop_order.append(pop_num)
                    pop_info[pop_num] = {'n_individuals': n_individuals, 'order': len(given_pop_order)}
            
            # Detectar início da seção de indivíduos
            if 'Inferred ancestry of individuals:' in line:
                reading_individuals = True
                next(f)  # Pular linha de cabeçalho
                continue
            
            # Ler dados dos indivíduos
            if reading_individuals:
                line = line.strip()
                if not line or line.startswith('--'):
                    break
                
                # Parse da linha: ID Label (%Miss) Pop: cluster1 cluster2 ...
                # Exemplo: "  1  HG00096    (0)    6 :  0.024 0.002 0.007 0.967"
                parts = line.split()
                if len(parts) < 4:
                    continue
                
                try:
                    ind_id = parts[0]
                    label = parts[1]
                    
                    # Encontrar o índice onde está o ":"
                    colon_idx = None
                    for i, part in enumerate(parts):
                        if ':' in part:
                            colon_idx = i
                            break
                    
                    if colon_idx is None:
                        continue
                    
                    # A população é o número ANTES do ":" (último número antes dele)
                    # Se o ":" está sozinho, pop é parts[colon_idx-1]
                    # Se o ":" está grudado (ex: "6:"), pop é parts[colon_idx] sem ":"
                    if parts[colon_idx] == ':':
                        pop = parts[colon_idx - 1]
                        cluster_start = colon_idx + 1
                    else:
                        pop = parts[colon_idx].rstrip(':')
                        cluster_start = colon_idx + 1
                    
                    # Os valores dos clusters vêm depois do ":"
                    cluster_values = [float(x) for x in parts[cluster_start:cluster_start+num_clusters]]
                    
                    if len(cluster_values) == num_clusters:
                        individuals.append({
                            'id': ind_id,
                            'label': label,
                            'pop': pop,
                            'clusters': cluster_values
                        })
                        populations.append(pop)
                        clusters.append(cluster_values)
                except (ValueError, IndexError) as e:
                    continue
    
    return individuals, num_clusters, given_pop_order, pop_info

def parse_pop_names(filename):
    """
    Lê arquivo com nomes das populações (uma por linha)
    Retorna dicionário mapeando número da população para nome
    """
    pop_names = {}
    try:
        with open(filename, 'r') as f:
            for i, line in enumerate(f, start=1):
                name = line.strip()
                if name:  # Ignora linhas vazias
                    pop_names[str(i)] = name
        print(f"Nomes de populações carregados: {len(pop_names)} populações")
        for pop_num, name in pop_names.items():
            print(f"  Pop {pop_num}: {name}")
    except FileNotFoundError:
        print(f"Erro: arquivo '{filename}' não encontrado")
        return None
    except Exception as e:
        print(f"Erro ao ler arquivo de nomes: {e}")
        return None
    
    return pop_names

def parse_pop_order(filename, pop_names=None):
    """
    Lê arquivo com ordem desejada das populações (uma por linha)
    Retorna lista com a ordem das populações
    Se pop_names for fornecido, valida os nomes
    """
    pop_order = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                name = line.strip()
                if name:  # Ignora linhas vazias
                    pop_order.append(name)
        print(f"\nOrdem customizada de populações carregada: {len(pop_order)} populações")
        print(f"Ordem: {' > '.join(pop_order)}")
        
        # Se temos nomes de populações, validar que os nomes na ordem existem
        if pop_names:
            valid_names = set(pop_names.values())
            for name in pop_order:
                if name not in valid_names:
                    print(f"Aviso: '{name}' na ordem não corresponde a nenhum nome de população")
    except FileNotFoundError:
        print(f"Erro: arquivo '{filename}' não encontrado")
        return None
    except Exception as e:
        print(f"Erro ao ler arquivo de ordem: {e}")
        return None
    
    return pop_order

def parse_colors(color_file_path):
    """
    Lê cores de um arquivo TXT.
    Cada linha pode conter uma cor em algum dos seguintes formatos:
      - Nome da cor (ex: red)
      - Código HEX (#56B4E9 ou 56B4E9)
      - RGB (ex: 230,159,0 ou (230,159,0))
    Linhas em branco ou comentários (linha começando com '# ') são ignoradas.
    """
    colors = []

    try:
        with open(color_file_path, 'r', encoding='utf-8-sig') as f:
            for line in f:
                color = line.strip()

                # Ignorar linhas vazias ou comentários explícitos
                if not color or color.startswith('# '):
                    continue

                try:
                    # HEX com ou sem '#'
                    if color.startswith('#') or all(c in '0123456789ABCDEFabcdef' for c in color):
                        if not color.startswith('#'):
                            color = '#' + color
                        colors.append(color)

                    # RGB com parênteses
                    elif color.startswith('(') and color.endswith(')'):
                        rgb = tuple(float(x.strip())/255 for x in color[1:-1].split(','))
                        colors.append(rgb)

                    # RGB simples (r,g,b)
                    elif ',' in color:
                        rgb = tuple(float(x.strip())/255 for x in color.split(','))
                        colors.append(rgb)

                    # Nome da cor
                    else:
                        colors.append(color)

                except Exception as e:
                    print(f"Aviso: cor inválida '{color}' ({e}), ignorada")

    except FileNotFoundError:
        print(f"Erro: arquivo de cores '{color_file_path}' não encontrado.")
        return None

    if not colors:
        print("Aviso: nenhum valor de cor válido encontrado. Usando paleta padrão.")
        return None

    return colors

def determine_cluster_color_mapping(individuals, num_clusters, custom_pop_order, pop_names):
    """
    Determina o mapeamento entre clusters e cores baseado na ancestralidade 
    predominante das primeiras populações da lista ordenada.
    
    Retorna um dicionário: {cluster_idx: color_idx}
    onde color_idx indica qual cor (na ordem do arquivo de cores) deve ser usada
    """
    if not custom_pop_order or not pop_names:
        # Se não há ordem customizada, usar mapeamento direto
        return {i: i for i in range(num_clusters)}
    
    # Criar mapeamento de nome para número de população
    name_to_num = {name: num for num, name in pop_names.items()}
    
    # Para cada população na ordem customizada, calcular a ancestralidade média por cluster
    cluster_mapping = {}
    assigned_clusters = set()
    
    print("\n=== Determinando mapeamento de cores ===")
    
    for color_idx, pop_name in enumerate(custom_pop_order):
        if color_idx >= num_clusters:
            break  # Já temos cores para todos os clusters
            
        pop_num = name_to_num.get(pop_name)
        if not pop_num:
            print(f"Aviso: população '{pop_name}' não encontrada")
            continue
        
        # Filtrar indivíduos desta população
        pop_individuals = [ind for ind in individuals if ind['pop'] == pop_num]
        
        if not pop_individuals:
            print(f"Aviso: nenhum indivíduo encontrado para população '{pop_name}'")
            continue
        
        # Calcular ancestralidade média por cluster para esta população
        cluster_means = [0.0] * num_clusters
        for ind in pop_individuals:
            for cluster_idx, value in enumerate(ind['clusters']):
                cluster_means[cluster_idx] += value
        
        cluster_means = [mean / len(pop_individuals) for mean in cluster_means]
        
        # Encontrar o cluster com maior ancestralidade média que ainda não foi atribuído
        sorted_clusters = sorted(range(num_clusters), key=lambda x: cluster_means[x], reverse=True)
        
        for cluster_idx in sorted_clusters:
            if cluster_idx not in assigned_clusters:
                cluster_mapping[cluster_idx] = color_idx
                assigned_clusters.add(cluster_idx)
                print(f"  {pop_name} -> Cluster {cluster_idx+1} (ancestralidade: {cluster_means[cluster_idx]:.3f}) -> Cor {color_idx+1}")
                break
    
    # Atribuir cores restantes aos clusters não mapeados
    unmapped_clusters = [i for i in range(num_clusters) if i not in cluster_mapping]
    for i, cluster_idx in enumerate(unmapped_clusters):
        color_idx = len(cluster_mapping) + i
        cluster_mapping[cluster_idx] = color_idx
        print(f"  Cluster {cluster_idx+1} não mapeado -> Cor {color_idx+1}")
    
    print("===================================\n")
    
    return cluster_mapping

def plot_structure_barplot(individuals, num_clusters, given_pop_order, pop_info,
                          output_file='structure_barplot.png', 
                          dpi=300, figsize=(20, 6), custom_colors=None,
                          pop_names=None, custom_pop_order=None):
    """
    Gera o gráfico de barras estilo CLUMPAK com amostras agrupadas por população.
    Agora com cores consistentes baseadas na ancestralidade predominante das primeiras populações.
    """

    # Determinar o mapeamento cluster -> cor
    cluster_to_color = determine_cluster_color_mapping(
        individuals, num_clusters, custom_pop_order, pop_names
    )

    # Determinar a ordem de exibição das populações
    display_order = []
    
    if custom_pop_order and pop_names:
        name_to_num = {name: num for num, name in pop_names.items()}
        for pop_name in custom_pop_order:
            if pop_name in name_to_num:
                display_order.append(name_to_num[pop_name])
            else:
                print(f"Aviso: '{pop_name}' não encontrado nos nomes de populações")
        for pop_num in given_pop_order:
            if pop_num not in display_order:
                display_order.append(pop_num)
                pop_label = pop_names.get(pop_num, f"Pop {pop_num}") if pop_names else f"Pop {pop_num}"
                print(f"Aviso: {pop_label} não estava na ordem customizada, adicionado ao final")
    else:
        display_order = given_pop_order

    # Ordenar indivíduos por população
    if display_order:
        pop_order_dict = {pop: idx for idx, pop in enumerate(display_order)}
        def sort_key(x):
            pop = x['pop']
            if pop in pop_order_dict:
                try:
                    ind_num = int(x['id'])
                except (ValueError, TypeError):
                    ind_num = 0
                return (pop_order_dict[pop], ind_num)
            else:
                try:
                    ind_num = int(x['id'])
                except (ValueError, TypeError):
                    ind_num = 0
                return (999, ind_num)
        individuals = sorted(individuals, key=sort_key)
    else:
        def sort_key(x):
            try:
                pop_num = int(x['pop']) if x['pop'] else 999
            except ValueError:
                pop_num = 999
            try:
                ind_num = int(x['id'])
            except (ValueError, TypeError):
                ind_num = 0
            return (pop_num, ind_num)
        individuals = sorted(individuals, key=sort_key)
    
    n_individuals = len(individuals)
    
    # Definir paleta de cores base
    if custom_colors and len(custom_colors) >= num_clusters:
        base_colors = custom_colors
    else:
        base_colors = list(plt.cm.Set3(np.linspace(0, 1, num_clusters)))
    
    # Aplicar o mapeamento de cores aos clusters
    colors = [None] * num_clusters
    for cluster_idx in range(num_clusters):
        color_idx = cluster_to_color.get(cluster_idx, cluster_idx)
        if color_idx < len(base_colors):
            colors[cluster_idx] = base_colors[color_idx]
        else:
            colors[cluster_idx] = base_colors[cluster_idx % len(base_colors)]
    
    # Criar figura
    fig, ax = plt.subplots(figsize=figsize)
    bottoms = np.zeros(n_individuals)
    
    for cluster_idx in range(num_clusters):
        values = [ind['clusters'][cluster_idx] for ind in individuals]
        ax.bar(range(n_individuals), values, bottom=bottoms, 
               color=colors[cluster_idx], width=1.0, edgecolor='none',
               label=f'Cluster {cluster_idx + 1}')
        bottoms += values
    
    # Linhas separadoras entre populações
    prev_pop = individuals[0]['pop']
    for i, ind in enumerate(individuals[1:], 1):
        if ind['pop'] != prev_pop:
            ax.axvline(x=i-0.5, color='black', linewidth=1, zorder=10)
            prev_pop = ind['pop']
    
    # Calcular posição média das populações
    pop_positions = {}
    for i, ind in enumerate(individuals):
        pop = ind['pop']
        if pop not in pop_positions:
            pop_positions[pop] = []
        pop_positions[pop].append(i)
    
    # Adicionar labels com ajuste automático (repelimento)
    text_positions = []
    if display_order:
        pops_iter = display_order
    else:
        pops_iter = sorted(pop_positions.keys(), key=lambda x: int(x) if str(x).isdigit() else x)
    
    for pop in pops_iter:
        if pop in pop_positions:
            positions = pop_positions[pop]
            mid_pos = (positions[0] + positions[-1]) / 2
            n_inds = len(positions)
            label_text = f"{pop_names.get(pop, f'Pop {pop}') if pop_names else f'Pop {pop}'}"
            
            # Ajustar posição vertical automaticamente se estiver muito próximo do anterior
            y_pos = -0.15
            for prev_x, prev_y in text_positions:
                if abs(mid_pos - prev_x) < 40:  # distância mínima em pixels
                    y_pos -= 0.05  # desloca um pouco para baixo
            
            text = ax.text(mid_pos, y_pos, label_text, 
                           ha='center', va='top', fontsize=10, fontweight='bold',
                           transform=ax.get_xaxis_transform())
            text_positions.append((mid_pos, y_pos))
    
    # Configurações do gráfico
    ax.set_xlim(-0.5, n_individuals - 0.5)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Ancestry proportion', fontsize=12, fontweight='bold')
    ax.set_xlabel('Individuals', fontsize=12, fontweight='bold')
    ax.set_title(f'STRUCTURE Analysis (K={num_clusters})', 
                fontsize=14, fontweight='bold', pad=20)
    
    ax.set_xticks([])
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), 
             frameon=True, fontsize=10)
    
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    print(f"Gráfico salvo em: {output_file}")
    
    return fig, ax

def main():
    parser = argparse.ArgumentParser(
        description='Gera gráfico de barras a partir do output do STRUCTURE com cores consistentes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplo de uso:
  python structure_barplot.py -i structure_output.txt -o meu_grafico.png
  python structure_barplot.py -i structure_output.txt -o output.png --dpi 600
  python structure_barplot.py -i structure_output.txt --figsize 30 8
  python structure_barplot.py -i structure_output.txt -c cores.txt
  python structure_barplot.py -i structure_output.txt -p pop_names.txt -l pop_order.txt -c cores.txt

IMPORTANTE: Para cores consistentes entre diferentes runs do STRUCTURE:
  - Use -p para fornecer os nomes das populações
  - Use -l para fornecer a ordem das populações (determina quais cores usar)
  - Use -c para fornecer o arquivo de cores
  As primeiras populações em -l receberão as primeiras cores em -c baseado na 
  ancestralidade predominante de cada população.
        """
    )
    
    parser.add_argument('-i', '--input', 
                       required=True,
                       dest='input_file',
                       help='Arquivo de output do STRUCTURE')
    parser.add_argument('-o', '--output', 
                       default='structure_barplot.png',
                       help='Nome do arquivo de saída (default: structure_barplot.png)')
    parser.add_argument('-p', '--popnames',
                       default=None,
                       help='Arquivo txt com nomes das populações (uma por linha, na ordem numérica)')
    parser.add_argument('-l', '--poporder',
                       default=None,
                       help='Arquivo txt com ordem customizada das populações (uma por linha, usando os nomes do -p)')
    parser.add_argument('-c', '--colors',
                       default=None,
                       help='Arquivo txt com cores customizadas (uma cor por linha). '
                            'Aceita nomes (red), HEX (#FF5733) ou RGB (255,87,51)')
    parser.add_argument('--dpi', 
                       type=int, 
                       default=300,
                       help='Resolução da imagem (default: 300)')
    parser.add_argument('--figsize', 
                       nargs=2, 
                       type=float, 
                       default=[20, 6],
                       help='Tamanho da figura em polegadas: largura altura (default: 20 6)')
    
    args = parser.parse_args()
    
    # Parse dos nomes das populações
    pop_names = None
    if args.popnames:
        pop_names = parse_pop_names(args.popnames)
        if pop_names is None:
            print("Continuando sem nomes customizados de populações...")
    
    # Parse da ordem das populações
    custom_pop_order = None
    if args.poporder:
        custom_pop_order = parse_pop_order(args.poporder, pop_names)
        if custom_pop_order is None:
            print("Continuando sem ordem customizada de populações...")
        elif not pop_names:
            print("Aviso: ordem de populações (-l) requer nomes de populações (-p). Ignorando ordem.")
            custom_pop_order = None
    
    # Parse das cores customizadas
    custom_colors = None
    if args.colors:
        custom_colors = parse_colors(args.colors)
        if custom_colors:
            print(f"Usando {len(custom_colors)} cores customizadas do arquivo")
    
    # Parse do arquivo
    print(f"Lendo arquivo: {args.input_file}")
    individuals, num_clusters, given_pop_order, pop_info = parse_structure_output(args.input_file)
    
    print(f"Total de indivíduos: {len(individuals)}")
    print(f"Número de clusters (K): {num_clusters}")
    print(f"Populações encontradas (ordem): {', '.join(given_pop_order)}")
    
    # Validar número de cores
    if custom_colors and len(custom_colors) < num_clusters:
        print(f"Aviso: apenas {len(custom_colors)} cores fornecidas para {num_clusters} clusters.")
        print("Cores faltantes serão geradas automaticamente.")
    
    # Gerar gráfico
    plot_structure_barplot(
        individuals, 
        num_clusters,
        given_pop_order,
        pop_info,
        output_file=args.output,
        dpi=args.dpi,
        figsize=tuple(args.figsize),
        custom_colors=custom_colors,
        pop_names=pop_names,
        custom_pop_order=custom_pop_order
    )
    
    print("Concluído!")


if __name__ == "__main__":
    main()