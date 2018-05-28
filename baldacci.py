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
