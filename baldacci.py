# We solve the GENPATH problem
def GENPATH(Delta, gamma, h, capacity, N, quantities, distance, direction):
    P = {}
    T = {}
    for k in N + [h]:
        P[k] = []
        T[k] = []
    T[h] = [{'path':[h],'cost':0,'lower_bound':0, "load":0, 'end':h}]

    count_paths = 0
    while True:
        costs = {}
        for k in T.keys():
            if len(T[k])>0:
                costs[k] = T[k][0]['cost']
        if len(costs) == 0:
            break
        min_costs = min(costs, key = costs.get)

        p_star = T[min_costs].pop(0)
        if not min_costs in P.keys():
            P[min_costs] = []
        P[min_costs].append(p_star)
        count_paths += 1
        # If too many paths, stop
        if count_paths==Delta:
            break
        # If path violates capacity, go to the next one
        if p_star['load'] > capacity/2.0:
            continue
        for n in N:
            if not (n in p_star['path']):
                if direction == 'right':
                    new_p = {'path':p_star['path'] + [n], 'cost':p_star['cost'] + distance[n][p_star['end']],
                            'lower_bound':p_star['cost'] + distance[n][p_star['end']], 'load': p_star['load'] + quantities[n],
                            'end':n}
                elif direction == 'left':
                    new_p = {'path':p_star['path'] + [n], 'cost':p_star['cost'] + distance[p_star['end']][n],
                            'lower_bound':p_star['cost'] + distance[p_star['end']][n], 'load': p_star['load'] + quantities[n],
                            'end':n}
                # Check if the new path has a cost too high
                if new_p['lower_bound'] >= gamma:
                    continue
                # Check if the new path has a load too high
                if new_p['load'] > capacity:
                    continue
                # Check if this new path is dominated by any path in P
                dominated = False
                for p in P[n]:
                    if (p['end'] == new_p['end']) and (p['cost'] <= new_p['cost']) and (set(p['path']) == set(new_p['path'])):
                        dominated = True
                        break
                if dominated:
                    continue
                # Check if the path is dominated by any path in T
                insertion_index = 0
                for i,p in enumerate(T[n]):
                    if (p['end'] == new_p['end']) and (p['cost'] <= new_p['cost']) and (set(p['path']) == set(new_p['path'])):
                        dominated = True
                        break
                    if (p['cost'] > new_p['cost']):
                        break
                    insertion_index = i+1
                if dominated:
                    continue
                # Append the path
                T[n].insert(insertion_index, new_p)
                # Delete dominated elements
                j = insertion_index + 1
                while j<len(T[n]):
                    p = T[n][j]
                    if (p['end'] == new_p['end']) and (p['cost'] > new_p['cost']) and (set(p['path']) == set(new_p['path'])):
                        T[n].pop(j)
                    else:
                        j += 1
    return P

def GENROUTE(Delta, gamma, h, capacity, N, quantities, distance):
    P_l = GENPATH(Delta, gamma, h, capacity, N, quantities, distance, direction = 'left')
    P_r = GENPATH(Delta, gamma, h, capacity, N, quantities, distance, direction = 'right')

    T = {}
    R = {}
    added = {}
    for n in N:
        added[n] = set((-1,-1))
        if len(P_l[n])>1 and len(P_r[n])>1:
            T[n] = [[(0,0),P_l[n][0]['cost']+P_r[n][0]['cost']]]
            added[n].add((0,0))
        else:
            T[n] = []
        R[n] = []

    valid_v = [0,0,0,0]
    while True:

        # Calculate costs
        costs = {}
        for n in N:
            if len(T[n])>0:
                costs[n] = T[n][0][1]
        if len(costs) == 0:
            break
        min_costs_n = min(costs, key = costs.get)
        min_cost = costs[min_costs_n]
        indices = T[min_costs_n].pop(0)[0]
        path_l = P_l[min_costs_n][indices[0]]
        path_r = P_r[min_costs_n][indices[1]]
        if min_cost> gamma:
            break
        total_load = path_l['load'] + path_r['load'] - quantities[min_costs_n]
        valid = True
        if total_load > capacity:
            valid = False
            valid_v[0] = valid_v[0]+1

        elif (np.min([path_l['load'],path_r['load']]) < total_load/2.0 or
            np.max([path_l['load'],path_r['load']]) > total_load/2.0+quantities[min_costs_n]):
            valid = False
            valid_v[1] = valid_v[1]+1

        elif (set(path_l['path']).intersection(set(path_r['path'])) != set([h,min_costs_n])):
            valid = False
            valid_v[2] = valid_v[2]+1
        else:
            for n in N:
                for r in R[n]:
                    if set(r['path']) == set(path_l['path']+path_r['path']):
                        valid = False
                        valid_v[3] = valid_v[3] + 1
                        break
                if not valid:
                    break
        if valid:
            R[min_costs_n].append({'path':path_l['path'][0:(len(path_l['path'])-1)]+list(reversed(path_r['path'])),
                                  'cost':path_l['cost']+path_r['cost'],
                                  'load':total_load,
                                  'median':min_costs_n,
                                  'indices':indices})

        new_route_1 = (indices[0]+1,indices[1])
        new_route_2 = (indices[0],indices[1]+1)
        # If routes do not exist, transform them into the first route
        if (indices[0]+1 >= len (P_l[min_costs_n])):
            new_route_1 = (0,0)
        if (indices[1]+1 >= len (P_r[min_costs_n])):
            new_route_2 = (0,0)
        new_routes = [new_route_1,new_route_2]
        new_costs = [P_l[min_costs_n][new_routes[0][0]]['cost']+P_r[min_costs_n][new_routes[0][1]]['cost'],
                     P_l[min_costs_n][new_routes[1][0]]['cost']+P_r[min_costs_n][new_routes[1][1]]['cost']]
        min_cost = np.min(new_costs)
        max_cost = np.max(new_costs)
        min_route = new_routes[np.argmin(new_costs)]
        max_route = new_routes[(np.argmin(new_costs)+1)%2]
        insert_index = 0
        # Check if the route has been added previously
        if not min_route in added[min_costs_n]:
            for i in range(len(T[min_costs_n])):
                cost = T[min_costs_n][i][1]
                if min_cost<cost:
                    break
                insert_index+=1
            T[min_costs_n].insert(insert_index,[min_route,min_cost])
            insert_index +=1
            added[min_costs_n].add(min_route)
        # Check if the route has been added previously
        if not max_route in added[min_costs_n]:
            for i in range(insert_index, len(T[min_costs_n])):
                cost = T[min_costs_n][i][1]
                if max_cost<cost:
                    break
                insert_index+=1
            T[min_costs_n].insert(insert_index,[max_route,max_cost])
            added[min_costs_n].add(max_route)

    # Verify that the routes are not always empty
    for n in N:
        if len(R[n]) == 0:
            R[n] = [{'path':[h,n,h], 'cost':distance[h][n] + distance[n][h], 'load': quantities[n]}]
    return R
