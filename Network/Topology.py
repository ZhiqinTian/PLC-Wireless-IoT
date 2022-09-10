import os
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt



class Topology:
    def __init__(self):
        self.topo = nx.Graph()
        self.nodes = []
        self.bone_paths = {}
        self.pos = {}
        self.base_area=set()
        self.areas = {}

    def get_pos(self):
        return self.pos

    def show(self):
        # pos = nx.spring_layout(self.topo)
        plt.figure(0)
        pos = self.pos
        if not pos:
            pos = nx.kamada_kawai_layout(self.topo)
        nx.draw(self.topo, pos=pos, node_color='red',
                edge_color='black', with_labels=False, node_size=50)
        node_labels = nx.get_node_attributes(self.topo, 'id')
        nx.draw_networkx_labels(
            self.topo, pos, font_size=5, labels=node_labels)
        edge_labels = nx.get_edge_attributes(self.topo, 'lenth')
        for edge in edge_labels:
            edge_labels[edge] = round(edge_labels[edge], 2)
        nx.draw_networkx_edge_labels(
            self.topo, pos, font_size=5, edge_labels=edge_labels)
        plt.savefig(os.getcwd()+"\\output_data\\topo.svg", format="svg")

    def add_node(self, node, position, load=1e8):
        self.nodes.append(node)
        self.base_area.add(node.id)
        pos = np.array(position)
        self.pos[node.id] = pos
        self.topo.add_node(node.id, id=node.id,
                           position=pos, load=load)

    def add_nodes(self, nodes, positions):
        for i in range(len(nodes)):
            self.add_node(nodes[i], positions[i])

    def set_load(self, node, load):
        self.topo.nodes[node.id]['load'] = load

    def set_loads(self, nodes, loads):
        for i in range(len(nodes)):
            self.topo.nodes[nodes[i].id]['load'] = loads[i]

    def add_line(self, node1, node2, lenth=None):
        if lenth is None:
            pos1 = self.topo.nodes[node1.id]['position']
            pos2 = self.topo.nodes[node2.id]['position']
            lenth = np.sqrt(np.sum(np.square(pos1-pos2)))
        self.topo.add_edge(node1.id, node2.id, lenth=lenth)

    def set_node_pos(self, node, position):
        pos = np.array(position)
        self.pos[node.id] = pos
        self.topo.nodes[node.id]['position'] = pos

    def set_line_lenth(self, node1, node2, lenth):
        self.topo.edges[node1.id, node2.id]['lenth'] = lenth

    def set_area(self, area_id):
        if area_id not in self.areas:
            self.areas[area_id]=set()
    
    def add_to_area(self,node,area_id):
        self.areas[area_id].add(node.id)
        if node.id in self.base_area:
            self.base_area.remove(node.id)

    def get_nodes_dis(self, node1, node2):
        pos1 = self.topo.nodes[node1.id]['position']
        pos2 = self.topo.nodes[node2.id]['position']
        dis = np.sqrt(np.sum(np.square(pos1-pos2)))
        return dis

    class Propagation:
        LOS = 0
        NLOS = 1

    def get_propagation(self, tx_node, rx_node):
        if (tx_node.id in self.base_area) and (rx_node.id in self.base_area):
            return self.Propagation.LOS
        for area_id in self.areas:
            if (tx_node.id in self.areas[area_id]) and (rx_node.id in self.areas[area_id]):
                return self.Propagation.LOS
        return self.Propagation.NLOS

    def bone_paths_refresh(self):
        self.bone_paths = dict(nx.all_pairs_shortest_path(self.topo))

    def units_struct_generate(self, tx_node, rx_node):

        def branch_level_struct(topo, level, black_list, connect_node_num, connect_node, unit):
            if len(unit['level_cable_num']) != level+1:
                unit['level_cable_num'].append(0)
                unit['level_load_type'].append([])
                unit['level_cable_type'].append([])
                unit['level_len'].append([])
                unit['next_connect'].append([])
                unit['level_node_load'].append([])
            neighbors = topo.neighbors(connect_node)
            next_node_num = 0
            black_list.append(connect_node)
            for next_node in neighbors:
                if next_node not in black_list:
                    unit['level_cable_num'][level] += 1
                    unit['level_load_type'][level].append(-1)
                    unit['level_cable_type'][level].append(1)
                    unit['level_len'][level].append(
                        topo.edges[connect_node, next_node]['lenth'])
                    unit['next_connect'][level].append(
                        connect_node_num)
                    unit['level_node_load'][level].append(
                        topo.nodes[next_node]['load'])
                    next_node_num += 1
                    if topo.degree(next_node) > 1:
                        branch_level_struct(
                            topo, level+1, black_list, next_node_num, next_node, unit)

        tx = tx_node.id
        rx = rx_node.id
        units_struct = []
        bone_path = self.bone_paths[rx][tx]
        unit_num = len(bone_path)
        for n in range(unit_num):
            unit = {}
            black_list = []
            branch_root = bone_path[n]
            branch_existence = 0
            root_degree = self.topo.degree(branch_root)
            # bone structure
            if n == 0:
                unit['bone1_len'] = 0
                unit['bone2_len'] = self.topo.edges[branch_root,
                                                    bone_path[n+1]]['lenth']/2
                black_list.append(bone_path[n+1])
                if root_degree > 1:
                    branch_existence = 1
            elif n == unit_num-1:
                unit['bone1_len'] = self.topo.edges[branch_root,
                                                    bone_path[n-1]]['lenth']/2
                unit['bone2_len'] = 0
                black_list.append(bone_path[n-1])
                if root_degree > 1:
                    branch_existence = 1
            else:
                unit['bone1_len'] = self.topo.edges[branch_root,
                                                    bone_path[n-1]]['lenth']/2
                unit['bone2_len'] = self.topo.edges[branch_root,
                                                    bone_path[n+1]]['lenth']/2
                black_list.append(bone_path[n-1])
                black_list.append(bone_path[n+1])
                if root_degree > 2:
                    branch_existence = 1
            unit['bone1_type'] = 1
            unit['bone2_type'] = 1

            # branch structure
            unit['root_load'] = self.topo.nodes[branch_root]['load']
            unit['bflag'] = 1
            unit['level_cable_num'] = []
            unit['level_load_type'] = []
            unit['level_cable_type'] = []
            unit['level_len'] = []
            unit['next_connect'] = []
            unit['level_node_load'] = []
            if branch_existence:
                unit['bflag'] = 0
                branch_level_struct(self.topo, 0, black_list,
                                    0, branch_root, unit)
            units_struct.append(unit)
        return unit_num, units_struct

    def load_change(self):
        pass
