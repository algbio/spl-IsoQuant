############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict

logger = logging.getLogger('IsoQuant')


class IntronCollector:
    def __init__(self, gene_info, delta=0):
        self.gene_info = gene_info
        self.known_introns = set(self.gene_info.intron_profiles.features)
        self.delta = delta
        # all intron counts
        # clustered intron counts
        self.clustered_introns = defaultdict(int)
        # how introns were corrected after clustering
        self.intron_correction_map = {}
        self.discarded_introns = set()

    def collect_introns(self, read_assignments):
        all_introns = defaultdict(int)
        for assignment in read_assignments:
            if not assignment.corrected_introns:
                continue
            for intron in assignment.corrected_introns:
                all_introns[intron] += 1
        return all_introns

    def construct_similar_intron_map(self, all_introns):
        ordered_introns = sorted(all_introns.keys())
        similar_intron_map = defaultdict(list)
        # create a dict of similar introns
        for i, intron in enumerate(ordered_introns):
            j = i + 1
            while j < len(ordered_introns) and abs(ordered_introns[j][0] - intron[0]) <= self.delta:
                if abs(ordered_introns[j][1] - intron[1]) <= self.delta:
                    similar_intron_map[intron].append(ordered_introns[j])
                    similar_intron_map[ordered_introns[j]].append(intron)
                j += 1
        return similar_intron_map

    def cluster_introns(self, all_introns, min_count):
        similar_intron_map = self.construct_similar_intron_map(all_introns)
        introns_sorted_by_counts = sorted([(v, k) for k, v in all_introns.items()], reverse=True)
        for count, intron in introns_sorted_by_counts:
            if intron in self.known_introns:
                # known intron is always added as trustworthy
                self.clustered_introns[intron] = count
            elif intron in similar_intron_map:
                # intron has a similar intron
                similar_introns = []
                for similar_intron in similar_intron_map[intron]:
                    if similar_intron in self.clustered_introns:
                        # if similar intron was already added to the cluster with the higher count
                        similar_introns.append((count, similar_intron))

                if similar_introns:
                    # take the best one as a substitute
                    substitute_intron = sorted(similar_introns, reverse=True)[0]
                    self.clustered_introns[substitute_intron[1]] += substitute_intron[0]
                    self.intron_correction_map[intron] = substitute_intron[1]
                else:
                    # no introns were found
                    self.clustered_introns[intron] = count
            elif count < min_count:
                self.discarded_introns.add(intron)
                continue
            else:
                self.clustered_introns[intron] = count

    def process(self, read_assignments, min_count):
        all_introns = self.collect_introns(read_assignments)
        self.cluster_introns(all_introns, min_count)
        logger.debug(self.intron_correction_map)

    def add_substitute(self, original_intron, substitute_intron):
        self.clustered_introns[substitute_intron] += self.clustered_introns[original_intron]
        del self.clustered_introns[original_intron]
        self.intron_correction_map[original_intron] = substitute_intron

    def discard(self, intron):
        self.discarded_introns.add(intron)
        if intron in self.clustered_introns:
            del self.clustered_introns[intron]

    def simplify_correction_map(self):
        all_introns = self.intron_correction_map.keys()
        for intron in all_introns:
            subs = self.intron_correction_map[intron]
            if subs not in self.intron_correction_map:
                continue
            while subs in self.intron_correction_map:
                subs = self.intron_correction_map[subs]
            self.intron_correction_map[intron] = subs


class IntronGraph:
    def __init__(self, params, gene_info, read_assignments):
        self.params = params
        self.params.min_novel_intron_count = 2
        self.params.graph_clustering_ratio = 0.2
        self.params.graph_clustering_distance = 50
        self.params.min_novel_isolated_intron_abs = 5
        self.params.min_novel_isolated_intron_rel = 0.2
        self.gene_info = gene_info
        self.read_assignments = read_assignments
        self.incoming_edges = defaultdict(set)
        self.outgoing_edges = defaultdict(set)
        self.intron_collector = IntronCollector(gene_info, params.delta)
        logger.debug("Collecting introns for %s" % self.gene_info.gene_db_list[0].id)
        self.intron_collector.process(read_assignments, self.params.min_novel_intron_count)
        self.construct()
        self.print_graph()
        self.simplify()
        self.print_graph()

    def add_edge(self, v1, v2):
        if v1 in self.intron_collector.intron_correction_map:
            v1 = self.intron_collector.intron_correction_map[v1]
        if v2 in self.intron_collector.intron_correction_map:
            v2 = self.intron_collector.intron_correction_map[v2]
        self.outgoing_edges[v1].add(v2)
        self.incoming_edges[v2].add(v1)

    def is_isolated(self, v):
        return not self.outgoing_edges[v] and not self.incoming_edges[v]

    # merge vertex to its substitute, remove if isolated
    def collapse_vertex(self, to_collapse, substitute_vertex):
        self.outgoing_edges[substitute_vertex].update(self.outgoing_edges[to_collapse])
        for i in self.outgoing_edges[to_collapse]:
            self.incoming_edges[i].remove(to_collapse)
            self.incoming_edges[i].add(substitute_vertex)

        self.incoming_edges[substitute_vertex].update(self.incoming_edges[to_collapse])
        for i in self.incoming_edges[to_collapse]:
            self.outgoing_edges[i].remove(to_collapse)
            self.outgoing_edges[i].add(substitute_vertex)

        self.intron_collector.add_substitute(to_collapse, substitute_vertex)

    def construct(self):
        logger.debug("Constructing for %s" % self.gene_info.gene_db_list[0].id)
        for assignment in self.read_assignments:
            if any(intron in self.intron_collector.discarded_introns for intron in assignment.corrected_introns):
                continue

            for i in range(len(assignment.corrected_introns) - 1):
                intron1 = assignment.corrected_introns[i]
                intron2 = assignment.corrected_introns[i + 1]
                self.add_edge(intron1, intron2)

    def simplify(self):
        logger.debug("Simplifying graph")
        # check all outgoing edges
        to_remove = set()
        logger.debug("Removing outgoing tips and bulges")
        for current_intron in self.outgoing_edges.keys():
            out_introns = self.outgoing_edges[current_intron]
            substitute_dict = self.collapse_vertex_set(out_introns)
            for i in substitute_dict.keys():
                to_remove.add(i)
                self.collapse_vertex(i, substitute_dict[i])

        # check all incoming edges
        logger.debug("Removing incoming tips and bulges")
        for current_intron in self.incoming_edges.keys():
            inc_introns = self.incoming_edges[current_intron]
            substitute_dict = self.collapse_vertex_set(inc_introns)
            for i in substitute_dict.keys():
                to_remove.add(i)
                self.collapse_vertex(i, substitute_dict[i])

        # check all isolated vertices
        isolated = set()
        logger.debug("Collapsing isolated introns")
        for intron in self.intron_collector.clustered_introns.keys():
            if self.is_isolated(intron):
                isolated.add(intron)
        substitute_dict = self.collapse_vertex_set(isolated)
        for i in substitute_dict.keys():
            to_remove.add(i)
            self.collapse_vertex(i, substitute_dict[i])

        # remove low covered isolated vertices
        logger.debug("Removing isolated introns")
        max_count = max(self.intron_collector.clustered_introns.values())
        count_cutoff = max(self.params.min_novel_isolated_intron_abs,
                           max_count * self.params.min_novel_isolated_intron_abs)
        for intron in isolated:
            if intron not in self.intron_collector.clustered_introns or intron in self.intron_collector.known_introns:
                # already removed or known
                continue
            if self.intron_collector.clustered_introns[intron] < count_cutoff:
                to_remove.add(intron)
                self.intron_collector.discard(intron)

        for i in to_remove:
            del self.outgoing_edges[i]
            del self.incoming_edges[i]

        self.intron_collector.simplify_correction_map()

    def collapse_vertex_set(self, vertex_set):
        if len(vertex_set) <= 1:
            return {}

        approved_set = set()
        substitute_dict = {}
        # get counts
        vertex_counts = sorted([(self.intron_collector.clustered_introns[i], i) for i in vertex_set], reverse=True)
        for count, vertex in vertex_counts:
            similar_vertices = []
            for i in approved_set:
                start_dist = abs(i[0] - vertex[0])
                end_dist = abs(i[1] - vertex[1])
                if start_dist < self.params.graph_clustering_distance and \
                        end_dist < self.params.graph_clustering_distance and \
                        count < self.intron_collector.clustered_introns[i] * self.params.graph_clustering_ratio:
                    # found similar intron among approved
                    similar_vertices.append((start_dist + end_dist, i))

            if not similar_vertices:
                # no similar introns
                approved_set.add(vertex)
            else:
                # selecting the most similar intron among approved
                substitute_vertex = sorted(similar_vertices)[0][1]
                substitute_dict[vertex] = substitute_vertex

        return substitute_dict

    def print_graph(self):
        logger.debug("Printing graph")
        for intron in sorted(self.intron_collector.clustered_introns.keys()):
            logger.debug("Intron %s, count %d -> %s" % (str(intron), self.intron_collector.clustered_introns[intron],
                                                        ",".join([str(x) for x in self.outgoing_edges[intron]])))

        for intron in sorted(self.intron_collector.clustered_introns.keys()):
            logger.debug("Intron %s, count %d <- %s" % (str(intron), self.intron_collector.clustered_introns[intron],
                                                        ",".join([str(x) for x in self.incoming_edges[intron]])))
