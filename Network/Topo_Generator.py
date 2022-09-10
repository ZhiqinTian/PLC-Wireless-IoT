from Network.Node import Node
from Network.Topology import Topology
import numpy as np

def indoor_senerio(N_room = 20,N_room_node = 6,Level_room = 3,N_room_hor = 10,
N_room_ver = 2,L_room = 10,H_room = 15,H_room_spacing=5):
    N_node = N_room*N_room_node
    nodes = Node.creat_nodes(N_node)
    topo = Topology()
    for i in range(N_room_ver):
        for j in range(N_room_hor):
            area_id=j+i*N_room_hor
            topo.set_area(area_id)
            node_id=j*N_room_node+i*N_room_node*N_room_hor
            pos = [j*L_room+L_room/2, i*(H_room+H_room_spacing)]
            topo.add_node(nodes[node_id], pos)
            topo.add_to_area(nodes[node_id],area_id)
            N_left = N_room_node-1
            last_level_node=[node_id]
            for k in range(Level_room):
                if N_left:
                    if k==Level_room-1:
                        N_l = N_left
                    else:
                        N_l = np.random.randint(min(N_left,N_room_node/2))+1
                    level_node=[]
                    for l in range(N_l):
                        node_id = node_id+1
                        level_node.append(node_id)
                        node_id_last=np.random.choice(last_level_node,1)[0]
                        pos_x=np.random.uniform(j*L_room+(l+0.2)*L_room/N_l,j*L_room+(l+0.8)*L_room/N_l)
                        pos_y=np.random.uniform((k+0.3)*H_room/Level_room+i*(H_room+H_room_spacing),(k+0.8)*H_room/Level_room+i*(H_room+H_room_spacing))
                        topo.add_node(nodes[node_id],[pos_x,pos_y])
                        topo.add_to_area(nodes[node_id],area_id)
                        topo.add_line(nodes[node_id_last], nodes[node_id])
                        N_left-=1
                    last_level_node=level_node
            if j>=1:
                topo.add_line(nodes[(j-1)*N_room_node+i*N_room_node*N_room_hor], nodes[j*N_room_node+i*N_room_node*N_room_hor])
        if i>=1:
            topo.add_line(nodes[(i-1)*N_room_node*N_room_hor], nodes[i*N_room_node*N_room_hor])
    
    return topo