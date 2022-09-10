from Network import *

nodes = Node.creat_nodes(3)
topo = Topology()
pos = [[0, 0], [0, 100], [10, 10]]
topo.set_nodes(nodes, pos)
topo.add_line(nodes[0],nodes[1])
topo.add_line(nodes[1],nodes[2])
topo.bone_paths_refresh()
unit_num, units_struct = topo.units_struct_generate(nodes[0], nodes[2])
print(unit_num)
for n in range(unit_num):
    print('---------------unit', n)
    for k in units_struct[n].keys():
        print(k, units_struct[n][k])
topo.show()
print(topo.get_pos())
