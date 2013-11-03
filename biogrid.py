#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import argparse
import httplib
import urllib
import csv
import StringIO
import itertools
import pydot
import networkx as nx


#if interactant in genes:
    #fillcolor = 'yellow'
    #if interactant in second_fact_genes:
        #if second_fact_genes[interactant] < -1:
            #fillcolor='magenta'
        #elif second_fact_genes[interactant] > 1:
            #fillcolor='cyan'
#elif interactant in second_fact_genes:
    #if second_fact_genes[interactant] < -1:
        #fillcolor='red'
    #elif second_fact_genes[interactant] > 1:
        #fillcolor='green'
BASE_URL = 'webservice.thebiogrid.org'


class BiogridConnection(object):
    def __init__(self, params={}):
        self.headers = {
            'Content-type': 'application/x-www-form-urlencoded',
            'Accept': 'text/plain'
        }
        self.params = {
            'searchNames': True,
            'start': 0,
            'enableCaching': True,
            'accesskey': '1ebc479d99a15ab52f9e3e231ea43bf9'
        }
        self.update_params(params)
        self.conn = httplib.HTTPConnection(BASE_URL)
    
    def update_params(self, p):
        self.params.update(p)
        for k in p:
            if p[k] is None:
                del p[k]

    def query(self, genes, include_input=True):
        qgenes = '|'.join(genes.keys())
        self.update_params({'geneList': qgenes})
        self.conn.request(
            'POST', '/interactions',
            body=urllib.urlencode(self.params),
            headers=self.headers
        )
        res = self.conn.getresponse()
        data = csv.reader(StringIO.StringIO(res.read()), delimiter='\t')
        interactions = []
        output_set = set()
        for row in data:
            output_set.add(row[7].upper())
            output_set.add(row[7].upper())
            try:
                g1, g2 = sorted([genes[row[7].upper()], genes[row[8].upper()]])
                interactions.append((g1, g2))
            except KeyError:
                pass
        interactions.sort()
        grouped_interactions = [(key, len(list(group))) for key, group in itertools.groupby(interactions)]
        if include_input:
            for gname, gene in genes.items():
                if gname not in output_set:
                    output_set.add(gname)
                    grouped_interactions.append(((gene, None), 0))
        return grouped_interactions


class Gene(object):
    def __init__(self, name, tags=[], contexts=[], meta={}):
        self.name = name.upper()
        self.tags = tags
        self.contexts = contexts
        self.meta = meta

    @classmethod
    def from_line(cls, line):
        l = line.strip('\n').split('\t')
        name = l[0]
        tags = map(str.strip, l[1].strip().split(';'))
        contexts = [c for c in map(str.strip, l[2].strip().split(';')) if c]
        meta = {}
        pairs = l[3].strip().split(';') 
        for p in pairs:
            meta.update(p)
        return cls(name, tags, contexts, meta)

    def __str__(self):
        return '\t'.join([
            self.name,
            ';'.join(self.tags),
            ';'.join(self.contexts),
            ';'.join([':'.join((k, v)) for k, v in self.meta])
        ])

    def __cmp__(self, other):
        return cmp(self.name, other.name)
                #and cmp(self.tags, other.tags) \
                #and cmp(self.functions, other.tags) \
                #and cmp(self.)

    def toJSON(self):
        return {
            'name': self.name,
            'groups': self.tags,
            'contexts': self.contexts,
            'meta': self.meta
        }


def parse_contexts_relationships(f, wl):
    rsh = [l.strip().split('\t') for l in f]
    if wl:
        o = []
        for r in rsh:
            if len(set(r) & set(wl)) == 2:
                o.append(r)
        return o
    return rsh

def link_contexts(graph, contexts, style={}):
    nodes, node_names = (set(), set())
    edges, edge_names = (set(), set())
    ctxt_style_default = {'shape': 'doubleoctagon', 'fillcolor': 'yellow', 'style': 'filled' }
    for pair in contexts:
        for p in pair:
            if p not in node_names:
                node_names.add(p)
                nodes.add(pydot.Node(p, **ctxt_style_default))
        if (pair[0], pair[1]) not in edge_names:
            edge_names.add((pair[0], pair[1]))
            edges.add(pydot.Edge(pair[0], pair[1], style='dashed', color='red', penwidth='2'))
    for node in nodes:
        graph.add_node(node)
    for edge in edges:
        graph.add_edge(edge)
    return graph

def parse_genes(f):
    genes = {}
    for line in f:
        gene = Gene.from_line(line)
        genes[gene.name] = gene
    return genes

def get_graph_nodes(graph):
    return [n.get_name() for n in graph.get_nodes()]

def draw_genes(graph, interactions, style={}):
    edges, edge_names = (set(), set())
    nodes = {}
    tag_style_default = { 'fillcolor': 'white', 'style': 'filled' }
    ctxt_style_default = {'shape': 'doubleoctagon', 'fillcolor': 'yellow', 'style': 'filled' }
    for interaction, evidence_count in interactions:
        tag_style = tag_style_default.copy()
        ctxt_style = ctxt_style_default.copy()
        for interactant in [i for i in interaction if i]:
            for tag in interactant.tags:
                if tag in style:
                    tag_style.update(style[tag])
            if interactant.name not in get_graph_nodes(graph):
                node = pydot.Node(interactant.name, **tag_style)
                nodes[interactant.name] = node
                graph.add_node(node)
                for ctxt in interactant.contexts:
                    if ctxt in style:
                        ctxt_style.update(style[ctxt])
                    if ctxt not in get_graph_nodes(graph):
                        node = pydot.Node(ctxt, **ctxt_style)
                        nodes[ctxt] = node
                        graph.add_node(node)
                    edges.add(pydot.Edge(nodes[interactant.name], nodes[ctxt]))
        if interaction[1] is not None and \
                (interaction[0].name, interaction[1].name) not in edge_names:
            edge_names.add((interaction[0].name, interaction[1].name))
            edges.add(pydot.Edge(nodes[interaction[0].name], nodes[interaction[1].name]))
                                 #penwidth=str(evidence_count)))
    for edge in edges:
        graph.add_edge(edge)
    return (graph, nodes)
    

def main(args):

    parameters = {
        'includeInteractors': args.includeInteractors,
        'includeInteractorInteractions': args.includeInteractorInteractions,
        'taxId': args.taxID,
    }
    if len(args.evidences) != 0:
        parameters.update({
            'includeEvidence': True,
            'evidenceList': '|'.join(args.evidences)
        })
    conn = BiogridConnection(parameters)

    genes = parse_genes(args.infile)
    
    grouped_interactions = conn.query(genes)
    
    if not args.keep_uniq:
        grouped_interactions = [i for i in grouped_interactions if i[1] > 0 or len(i[0][0].contexts) > 0]

    print [i[0][0] for i in grouped_interactions]
    if args.tag:
        pass
        #grouped_interactions = [i for i in grouped_interactions
                #if i[0][0].tags ]

    if args.context_whitelist:
        grouped_interactions = [i for i in grouped_interactions if len(set(i[0][0].contexts) & set(args.context_whitelist)) >= 1]
        for g in grouped_interactions:
            for i in g[0]:
                if i is not None:
                    i.contexts = list(set(i.contexts) & set(args.context_whitelist))

    graph = pydot.Dot(graph_type='graph', outputorder='nodeslast', splines='true')
    #graph = pydot.Dot(graph_type='graph', outputorder='edgesfirst')
    #graph = pydot.Dot(graph_type='graph', overlap='scale')
    #graph = pydot.Dot(graph_type='graph', pack='true')

    tag_style = {
        'cort_only': {'fillcolor': 'blue'},
        't3_only': {'fillcolor': 'magenta'},
        'antagonism': {'fillcolor': 'red'},
        'potentiation': {'fillcolor': 'green'}
    }
    graph, nodes = draw_genes(graph, grouped_interactions, style=tag_style)

    contexts = None
    if args.contexts:
        contexts = parse_contexts_relationships(args.contexts, args.context_whitelist)
        graph = link_contexts(graph, contexts)

    graph.write(args.outfile+'.dot')
    #graph.write_png('090913_cort_network.png', prog='twopi')
    graph.write_png(args.outfile+'.png', prog='fdp')
    #graph.write_png('090913_cort_network.png', prog='neato')
    #graph.write_png('090913_cort_network.png', prog='dot')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--infile', dest='infile',
        type=argparse.FileType('r'),
        help='Input file'
    )
    parser.add_argument(
        '-o', '--outfile', dest='outfile',
        help='Output file'
    )
    parser.add_argument(
        '-F', '--contexts', dest='contexts',
        type=argparse.FileType('r'),
        help='Interactions between functions/contexts/pathologies'
    )
    parser.add_argument(
        '-U', '--includeInteractors', dest='includeInteractors',
        action='store_true',
        default=False,
        help='Include first order interactors'
    )
    parser.add_argument(
        '-I', '--includeInteractorInteractions', dest='includeInteractorInteractions',
        action='store_true',
        default=False,
        help='Include interactions between first order interactors'
    )
    parser.add_argument(
        '-T', '--tag', dest='tag',
        type=type,
        default=None,
        help='Keep only interactions involving genes labeled with the specified tag'
    )
    parser.add_argument(
        '-e', '--evidences', dest='evidences',
        nargs='*',
        default=[],
        #default=['Biochemical Activity', 'Dosage Growth Defect', 'Dosage Lethality',
            #'Dosage Rescue', 'Negative Genetic', 'Phenotypic Enhancement', 'Phenotypic Suppression',
            #'Positive Genetic', 'Synthetic Growth Defect', 'Synthetic Haploinsufficiency', 'Synthetic Lethality'],
        help='Evidence types to consider'
    )
    parser.add_argument(
        '-t', '--taxID', dest='taxID',
        default='All',
        help='Taxonomic ID (9606 = H. sapiens; 8364 = X. tropicalis)'
    )
    parser.add_argument(
        '-u', '--keep_uniq', dest='keep_uniq',
        action='store_true',
        default=False,
        help='Do not discard genes not interacting with anything'
    )
    parser.add_argument(
        '-w', '--context_whitelist', dest='context_whitelist',
        nargs='*',
        default=[],
        help='Discard genes not belonging to at least these contexts'
    )
    parser.add_argument(
        '-O', '--get_non_redundant_interactors', dest='get_non_redundant_interactors',
        action='store_true',
        default=False,
        help='help'
    )
    main(parser.parse_args())
