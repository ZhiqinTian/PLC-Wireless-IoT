from Sim.Simulator import *
from ..Net_Protocol import Net_Route
from ..Packet import *
from ..Address import Addresses
from Sim.Time import *
import copy


class LOADng:

    class State:
        INVALID = 0
        VALID = 1
        IN_SEARCH = 2

    class Error:
        NO_ROUTE = 0
        TX_ERROR = 1

    def __init__(self, params=None):
        self.node = None
        self.rreq_id = 0
        self.rreq_id_cache = []

        # parameters of LOADng protocol algorithm
        if not params:
            params = self.Helper.LOADng_params()
        self.life_time = params.life_time
        self.bad_link_life_time = params.bad_link_life_time
        self.expire_time = params.expire_time

        # routing table and request queue
        self.routing_table = self.Routing_Table(
            self.life_time, self.bad_link_life_time)
        self.request_queue = self.Request_Queue(self.expire_time)

    def route_output(self, np_header):
        source = np_header.source
        destination = np_header.destination
        route = Net_Route()
        route.source = source
        route.destination = destination
        route_entry = self.routing_table.lookup(destination)
        if route_entry:
            if route_entry.state == self.State.VALID:
                route.next_hop = route_entry.next_hop
                route.tx_params = route_entry.tx_params
                self.routing_table.update_life_time(
                    destination, self.life_time.sim_time)
                self.routing_table.update_life_time(
                    route.next_hop, self.life_time.sim_time)
                return route
        route.next_hop = Addresses.LOOPBACK
        return route

    def deferred_route_output(self, net_packet):
        queue_entry = self.Request_Queue.Entry(net_packet)
        self.request_queue.add_entry(queue_entry)
        self.routing_table.refresh()
        destination = net_packet.header.destination
        route_entry = self.routing_table.lookup(destination)
        if (not route_entry) or (route_entry and (route_entry.state != self.State.IN_SEARCH)):
            self.send_rreq(destination)

    def route_input(self, net_packet, fr, to, net_send_data, net_receive_data, net_forward):
        source = net_packet.header.source
        destination = net_packet.header.destination
        payload = net_packet.payload
        if fr == Addresses.LOOPBACK:
            self.deferred_route_output(net_packet)
            return
        if to == source:
            print('Own packet. Drop.')
            return
        self.routing_table.refresh()
        self.routing_table.update_life_time(
            source, self.life_time.sim_time)
        self.routing_table.update_life_time(
            fr, self.life_time.sim_time)
        if payload.type != 'data_packet':
            msg_new, dst_new = self.receive_process(
                net_packet, fr, to, net_send_data, net_forward)
            if msg_new:
                if destination == Addresses.BROADCAST:
                    if dst_new:
                        net_send_data(to, dst_new, msg_new)
                        return
                    net_send_data(to, destination)
                else:
                    header_new = copy.copy(net_packet.header)
                    net_packet_new = Net_Packet(header_new, msg_new)
                    self.forwarding(net_packet_new, net_forward)
            return
        if (to == destination) or (destination == Addresses.BROADCAST):
            net_receive_data(net_packet)
            return
        self.forwarding(net_packet, net_forward)

    def forwarding(self, net_packet, net_forward):
        source = net_packet.header.source
        destination = net_packet.header.destination
        route = Net_Route()
        route.source = source
        route.destination = destination
        self.routing_table.refresh()
        route_entry = self.routing_table.lookup(destination)
        if route_entry:
            if route_entry.state == self.State.VALID:
                route.next_hop = route_entry.next_hop
                route.tx_params = route_entry.tx_params
                self.routing_table.update_life_time(
                    destination, self.life_time.sim_time)
                self.routing_table.update_life_time(
                    route.next_hop, self.life_time.sim_time)
                forward_result = net_forward(route, net_packet)
                if not forward_result:
                    self.route_error(self.Error.TX_ERROR, source,
                                     destination, route.next_hop)
                return
        self.route_error(self.Error.NO_ROUTE, source, destination, None)

    def send_rreq(self, destination):
        rreq = self.RREQ(self.node.net_device.address, destination)
        self.routing_table.refresh()
        route_entry = self.routing_table.lookup(destination)
        if route_entry:
            route_entry.state = self.State.IN_SEARCH
        else:
            new_entry = self.routing_table.Entry()
            new_entry.destination = destination
            new_entry.state = self.State.IN_SEARCH
            self.routing_table.add_entry(new_entry)
        self.rreq_id = self.rreq_id+1
        src = self.node.net_device.address
        dst = Addresses.BROADCAST
        rreq.id = self.rreq_id
        Simulator.schedule_now(
            self.node.net_protocol.send_data, src, dst, rreq)

    def send_rrep(self, rreq):
        destination = rreq.origin
        rrep = self.RREP(rreq.origin, rreq.destination)
        src = self.node.net_device.address
        dst = destination
        Simulator.schedule_now(
            self.node.net_protocol.send_data, src, dst, rrep)

    def route_error(self, error_type, error_packet_source, error_packet_destination, next_hop):
        if error_type == self.Error.TX_ERROR:
            self.routing_table.turn_bad(next_hop)
            self.routing_table.turn_bad(error_packet_destination)
        self.send_rerr(error_type, error_packet_source,
                       error_packet_destination)

    def send_rerr(self, error_type, error_packet_source, error_packet_destination):
        src = self.node.device.address
        dst = error_packet_source
        rerr = self.RERR(src, dst)
        rerr.error_type = error_type
        rerr.unreachable_dst = error_packet_destination
        Simulator.schedule_now(
            self.node.net_protocol.send_data, src, dst, rerr)

    def send_packet_from_queue(self, destination):
        net_forward = self.node.net_protocol.forward
        self.request_queue.refresh()
        while True:
            net_packet = self.request_queue.dequeue(destination)
            if not net_packet:
                break
            self.forwarding(net_packet, net_forward)

    def receive_process(self, net_packet, fr, to):
        LOADng_msg = net_packet.payload
        if LOADng_msg.type == 'rreq':
            msg_new, dst_new = self.receive_rreq(LOADng_msg, fr, to)
        elif LOADng_msg.type == 'rrep':
            msg_new, dst_new = self.receive_rrep(LOADng_msg, fr, to)
        elif LOADng_msg.type == 'rerr':
            msg_new, dst_new = self.receive_rerr(LOADng_msg, fr, to)
        return msg_new, dst_new

    def receive_rreq(self, rreq, fr, to):
        origin = rreq.origin
        destination = rreq.destination
        id_ = rreq.id
        rreq_id = (origin, id_)
        if rreq_id in self.rreq_id_cache:
            return None, None
        else:
            self.rreq_id_cache.append(rreq_id)
        hops = rreq.hops+1

        # create reverse route to the neighbor where the rreq come from
        self.routing_table.refresh()
        route_neighbor = self.routing_table.lookup(fr)
        if not route_neighbor:
            route_neighbor = self.routing_table.Entry()
            self.routing_table.add_entry(route_neighbor)
        route_neighbor.destination = fr
        route_neighbor.next_hop = fr
        route_neighbor.hops = 1
        route_neighbor.life_time = Simulator.now().sim_time+self.life_time.sim_time
        route_neighbor.state = self.State.VALID

        # create reverse route to the origin where the rreq come from
        route_origin = self.routing_table.lookup(origin)
        if not route_origin:
            route_origin = self.routing_table.Entry()
            self.routing_table.add_entry(route_origin)
        route_origin.destination = origin
        route_origin.next_hop = fr
        route_origin.hops = hops
        route_origin.life_time = Simulator.now().sim_time+self.life_time.sim_time
        route_origin.state = self.State.VALID

        # send rrep if the destination of the rreq is this router
        if destination == to:
            self.send_rrep(rreq)
            return None, None

        # if it has an active route to the destination
        route_dst = self.routing_table.lookup(destination)
        if route_dst:
            if route_dst.state == self.State.VALID:
                dst_new = destination
        rreq_new = copy.copy(rreq)
        rreq_new.hops = hops
        return rreq_new, dst_new

    def receive_rrep(self, rrep, fr, to):
        origin = rrep.origin
        destination = rrep.destination
        hops = rrep.hops+1
        # create reverse route to the neighbor where the rreq come from
        self.routing_table.refresh()
        route_neighbor = self.routing_table.lookup(fr)
        if not route_neighbor:
            route_neighbor = self.routing_table.Entry()
            self.routing_table.add_entry(route_neighbor)
        route_neighbor.destination = fr
        route_neighbor.next_hop = fr
        route_neighbor.hops = 1
        route_neighbor.life_time = Simulator.now().sim_time+self.life_time.sim_time
        route_neighbor.state = self.State.VALID

        # create reverse route to the destination where the rrep come from
        route_dst = self.routing_table.lookup(destination)
        if not route_dst:
            route_dst = self.routing_table.Entry()
            self.routing_table.add_entry(route_dst)
        route_dst.destination = destination
        route_dst.next_hop = fr
        route_dst.hops = hops
        route_dst.life_time = Simulator.now().sim_time+self.life_time.sim_time
        route_dst.state = self.State.VALID

        # if the origin of the rrep's rreq is this router
        if origin == to:
            self.send_packet_from_queue(destination)
            return None, None

        rrep_new = copy.copy(rrep)
        rrep_new.hops = hops
        return rrep_new, origin

    def receive_rerr(self, rerr, fr, to):
        dst = rerr.destination
        unreachable_dst = rerr.unreachable_dst
        self.routing_table.refresh()
        route_entry = self.routing_table.lookup(unreachable_dst)

        if route_entry:
            if route_entry.next_hop == fr:
                self.routing_table.turn_bad(unreachable_dst)

        if dst == to:
            return None, None

        rerr_new = copy.copy(rerr)
        return rerr_new, dst

    class Helper:

        def __init__(self):
            self.LOADng_params = self.LOADng_params()

        def set(self, params):
            self.LOADng_param = params

        def create(self):
            return LOADng(self.LOADng_param)

        class LOADng_params:
            def __init__(self):
                self.life_time = Second(10)
                self.bad_link_life_time = Second(5)
                self.expire_time = Second(5)

    class Routing_Table:
        def __init__(self, life_time, bad_link_life_time):
            self.table = dict()
            self.life_time = life_time
            self.bad_link_life_time = bad_link_life_time

        def lookup(self, destination):
            if destination.value in self.table:
                return self.table[destination.value]
            return None

        def add_entry(self, table_entry):
            dst = table_entry.destination.value
            self.table[dst] = table_entry

        def delete_entry(self, destination):
            dst = destination.value
            if dst in self.table:
                self.table.pop(dst)

        def update_life_time(self, destination, delta_life_time):
            dst = destination.value
            if dst in self.table:
                table_entry = self.table[dst]
                if table_entry.state == LOADng.State.VALID:
                    table_entry.life_time = table_entry.life_time+delta_life_time

        def turn_bad(self, destination):
            dst = destination.value
            if self.table[dst].state == LOADng.State.VALID:
                self.table[dst].state = LOADng.State.INVALID
                self.table[dst].life_time = self.life_time.sim_time + \
                    Simulator.now().sim_time

        def refresh(self):
            for dst in self.table:
                if self.table[dst].life_time < Simulator.now().sim_time:
                    if self.table[dst].state == LOADng.State.INVALID:
                        self.table.pop(dst)
                    elif self.table[dst].state == LOADng.State.VALID:
                        self.turn_bad(dst)

        class Entry:
            def __init__(self):
                self.destination = None
                self.next_hop = None
                self.hops = None
                self.state = None
                self.life_time = None
                self.tx_params = None

    class Request_Queue:
        def __init__(self, expire_time=None):
            self.queue = []
            self.expire_time = expire_time

        def add_entry(self, entry):
            entry.expire_time = self.expire_time.sim_time + Simulator.now().sim_time
            self.queue.append(entry)

        def dequeue(self, destination):
            for i in range(len(self.queue)):
                if self.queue[i].net_packet.header.destination == destination:
                    net_packet = self.queue[i].net_packet
                    self.queue.pop(i)
                    return net_packet
            return None

        def refresh(self):
            while True:
                if self.queue[0].expire_time < Simulator.now().sim_time:
                    self.queue.pop(0)
                else:
                    break

        class Entry:
            def __init__(self, net_packet=None):
                self.net_packet = net_packet
                self.expire_time = None

    class RREQ:
        def __init__(self, origin, destination):
            self.type = 'rreq'
            self.id = None
            self.origin = origin
            self.destination = destination
            self.route_cost = None
            self.hops = 0
            self.hop_limit = None

    class RREP:
        def __init__(self, origin, destination):
            self.type = 'rrep'
            self.origin = origin
            self.destination = destination
            self.hops = 0

    class RREP_ACK:
        def __init__(self,):
            self.type = 'rrep_ack'

    class RERR:
        def __init__(self, source, destination):
            self.type = 'rerr'
            self.error_type = None
            self.unreachable_dst
            self.hops_limit = None
            self.source = source
            self.destination = destination
