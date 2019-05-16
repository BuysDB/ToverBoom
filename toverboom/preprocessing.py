#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import networkx as nx

"""
Convert from '(1,2)' to (1,2) keys
"""
def convertFromStringKeys(graph):
    mapping = {}
    for key in graph:
        clone,tp = key.replace('(','').replace(')','').split(',')
        clone = int(clone)
        tp = int(tp)
        mapping[key] = (clone, tp)
    return nx.relabel_nodes(graph,mapping)
