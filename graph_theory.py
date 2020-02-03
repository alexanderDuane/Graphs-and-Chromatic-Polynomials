

class Graph:

    chromatic_polynomial = None

    def __init__(self, adjacency, cache=None,level=1):
        self.adjacency = dict(adjacency)
        self.level = level
        self.cache = cache

    @property
    def vertices(self):
        '''
        Set return
        '''
        return set(self.adjacency.keys())
    
    @property
    def edge_count(self):
        '''
        Integer return

        '''
        edge_count = 0
        for x in self.adjacency.values():
            edge_count += len(x)
        return edge_count
    
    @property
    def is_tree(self):
        '''
        Boolean return
        '''
        if (self.edge_count/2) == len(self.vertices) - 1:
            return True
        else:
            return False

    @property
    def is_connected(self):
        '''
        Boolean return
        '''
        vertex = tuple(self.vertices)[0]
        if self.traverse({vertex},vertex) == self.vertices:
            return True
        else:
            return False

    @property
    def best_nodes_for_deletion_contraction(self):
        '''
        Returns the two connected nodes of the highest degree possible where deleting the edge between them does not result in disconnecting the graph.
        '''

        decending_degree_list = sorted([v for v in self.vertices if len(self.adjacency[v]) > 1], key= lambda x: len(self.adjacency[x]), reverse=True)        
       
        if decending_degree_list == []:
            return False, False
        else:
            for vertex1 in decending_degree_list:

                vertex2_candidates = sorted([v for v in self.adjacency[vertex1] if len(self.adjacency[v]) > 1], key= lambda x: len(self.adjacency[x]), reverse=True)        
                    
                for vertex2 in vertex2_candidates:
                    
                    deletion_graph = self.deletion_graph(vertex1,vertex2)
                    if deletion_graph.is_connected:

                        contraction_graph = self.contraction_graph(vertex1,vertex2)

                        return vertex1, vertex2
 
    @property
    def deletion_contraction(self):
        '''
        Function return.
        Returns chromatic polynomial using the deletion-contraction method

        '''
        if self.is_connected == False:
            print('does not exist: graph is not connected')
            return None
    
        if self.is_tree:
            
            return lambda n: n*(n-1)**(len(self.vertices) - 1)
        
        else:

            vertex1, vertex2 = self.best_nodes_for_deletion_contraction
            deletion_poly    = self.deletion_graph(vertex1, vertex2).deletion_contraction
            contraction_poly = self.contraction_graph(vertex1, vertex2).deletion_contraction      
            return lambda G: deletion_poly(G) - contraction_poly(G)

    @property
    def chromatic_number(self):
        '''
        Int return
        The minimum number of colors needed for a coloring of self
        '''
        if self.chromatic_polynomial == None:
            self.glue()

        chromatic_number = 1

        while self.chromatic_polynomial(chromatic_number) <= 0:
            chromatic_number += 1
        return chromatic_number

    @staticmethod
    def subgraph_of_common_vertex(graph, vertex):
        '''
        Graph return
        Returns the subgraph of 'graph' where every vertex is connected to 'vertex'.
        '''
        v_subset = graph.adjacency[vertex].union({vertex})
        
        return Graph({v : {x for x in graph.adjacency[v] if x in v_subset} for v in v_subset})

    @staticmethod
    def subgraph(graph,vertices):
        '''
        Graph return
        Returns the subgraph of 'graph' where every vertex is in the set 'vertices'.
        '''
        return Graph({v : {x for x in graph.adjacency[v] if x in vertices} for v in vertices})
        
    @staticmethod
    def order_list_of_children(graph,vertex):
        '''
        List return
        Returns list of the set order of children for all vertices of 'vertex' in 'graph'
        '''
        return [len(graph.adjacency[children]) for children in graph.adjacency[vertex]]
        
    @staticmethod
    def degree_count_dictionary(graph):
        '''
        Dictionary return
        Returns a dictionary of the number of vertices of a degree paired to the degree number itself
        '''
        degree_count_dictionary = {}

        for vertex, connected_nodes in graph.adjacency.items():
            if len(connected_nodes) in degree_count_dictionary.keys():
                degree_count_dictionary[len(connected_nodes)].add(vertex)
            else:
                degree_count_dictionary[len(connected_nodes)] = {vertex}
        return degree_count_dictionary

    def traverse(self, vertex_set, vertex):
        connected = self.adjacency[vertex]
        to_be_traversed = set([v for v in connected if v not in vertex_set])
        
        if len(to_be_traversed) == 0:
            return vertex_set
        else:
            for v in to_be_traversed:
                vertex_set = self.traverse(vertex_set.union(to_be_traversed),v)
            return vertex_set
      
    def deletion_graph(self,vertex1,vertex2):
        
        adjacency = dict(self.adjacency)
        adjacency[vertex2] = set([v for v in adjacency[vertex2] if v != vertex1])
        adjacency[vertex1] = set([v for v in adjacency[vertex1] if v != vertex2])

        return Graph(adjacency=adjacency,level=self.level + 1)
    
    def contraction_graph(self,vertex1,vertex2):

        adjacency = dict(self.adjacency)
        adjacency[vertex2] = set([v for v in adjacency[vertex2] if v != vertex1])
        adjacency[vertex1] = set([v for v in adjacency[vertex1] if v != vertex2]).union(adjacency[vertex2])

        for node in adjacency[vertex2]:
            adjacency[node] = set([v for v in adjacency[node] if v != vertex2]).union({vertex1})
            
        adjacency[vertex1].union(adjacency[vertex2])
        del adjacency[vertex2]

        return Graph(adjacency=adjacency,level=self.level + 1)
    
    def glue(self):
        vertices_remaining = [v for v in self.vertices]
        vertex_1 = vertices_remaining.pop()

        initial_graph = self.subgraph_of_common_vertex(self,vertex_1)
        initial_graph.chromatic_polynomial = initial_graph.deletion_contraction

        self.chromatic_polynomial = lambda G: initial_graph.deletion_contraction(G)*self._glue(initial_graph, vertices_remaining)(G)

    def _glue(self, graph_main, vertices_remaining):

        if vertices_remaining:

            vertex_2 = [v for v in graph_main.vertices if v in vertices_remaining][0]
            vertices_remaining = [v for v in vertices_remaining if v != vertex_2]

            subgraph_vertex_2       = self.subgraph_of_common_vertex(self,vertex_2)
            intersection_subgraph   = self.subgraph(self,graph_main.vertices.intersection(subgraph_vertex_2.vertices))
                
            subgraph_vertex_2.chromatic_polynomial      = subgraph_vertex_2.deletion_contraction
            intersection_subgraph.chromatic_polynomial  = lambda G: max(intersection_subgraph.deletion_contraction(G),1)

            graph_main = self.subgraph(self,graph_main.vertices.union(subgraph_vertex_2.vertices))

            #alpha = self._glue(graph_main, vertices_remaining)


            return lambda G: self._glue(graph_main, vertices_remaining)(G)*subgraph_vertex_2.chromatic_polynomial(G)/intersection_subgraph.chromatic_polynomial(G)
            #return lambda G: alpha(G)*subgraph_vertex_2.chromatic_polynomial(G)/intersection_subgraph.chromatic_polynomial(G)
        
        else:
            
            return lambda G: 1

    def are_isomorphic(self,graph_1,graph_2):

        degree_count_dict1, degree_count_dict2 = (self.degree_count_dictionary(graph=G) for G in (graph_1,graph_2))

        if degree_count_dict1.keys() != degree_count_dict2.keys():
            return False

        if {key: len(values) for key,values in degree_count_dict1.items()} != {key: len(values) for key,values in degree_count_dict2.items()}:
            return False
        
        highest_degree = sorted(degree_count_dict1.keys(), reverse=True)[0]
        
        vertex_1 = degree_count_dict1[highest_degree].pop()
        
        vertex_2_candidates = [vertex for vertex in degree_count_dict2[highest_degree] if 
                               set(self.order_list_of_children(graph_2,vertex)) == set(self.order_list_of_children(graph_1,vertex_1))]

        for vertex_2 in vertex_2_candidates:
            is_isomorphic = self._are_isomorphic(graph_1,graph_2,dict(degree_count_dict1),dict(degree_count_dict2),vertex_1,vertex_2)
            if is_isomorphic:
                return True


        return False

    def _are_isomorphic(self, graph_1,graph_2,degree_count_dict1_orig,degree_count_dict2_orig,vertex_1,vertex_2,isomorphism_orig={}):
        isomorphism = dict(isomorphism_orig)

        for v in [v for v in graph_1.adjacency[vertex_1] if v in isomorphism.keys()]:
            if isomorphism[v] not in graph_2.adjacency[vertex_2]:
                return False
                                     
        isomorphism[vertex_1] = vertex_2
        
        degree_count_dict1 = dict(degree_count_dict1_orig)
        degree_count_dict2 = dict(degree_count_dict2_orig)
        
        if degree_count_dict1:

            highest_degree = sorted(degree_count_dict1.keys(), reverse=True)[0]
            vertex_11 = degree_count_dict1[highest_degree].pop()
            
            vertex_22_candidates = [vertex for vertex in degree_count_dict2[highest_degree] if 
                                   (set(self.order_list_of_children(graph_2,vertex)) == set(self.order_list_of_children(graph_1,vertex_11))) and (vertex not in isomorphism.values())]
            
            for vertex_22 in vertex_22_candidates:
                is_isomorphic = self._are_isomorphic(graph_1,graph_2,dict(degree_count_dict1),dict(degree_count_dict2),vertex_11,vertex_22,isomorphism)
                if is_isomorphic:
                    return True

            return False
        else:

            return True

class SampleGraphs:            
    america = {
                'Alabama':{'Florida','Georgia','Mississippi','Tennessee'},
                'Arizona':{'California','Colorado','Nevada','New_Mexico','Utah'},
                'Arkansas':{'Louisiana','Mississippi','Missouri','Oklahoma','Tennessee','Texas'},
                'California':{'Arizona','Nevada','Oregon'},
                'Colorado':{'Arizona','Kansas','Nebraska','New_Mexico','Oklahoma','Utah','Wyoming'},
                'Connecticut':{'Massachusetts','New_York','Rhode_Island'},
                'Delaware':{'Maryland','New_Jersey','Pennsylvania'},
                'Florida':{'Alabama','Georgia'},
                'Georgia':{'Alabama','Florida','North_Carolina','South_Carolina','Tennessee'},
                'Idaho':{'Montana','Nevada','Oregon','Utah','Washington','Wyoming'},
                'Illinois':{'Indiana','Iowa','Kentucky','Missouri','Wisconsin'},
                'Indiana':{'Illinois','Kentucky','Michigan','Ohio'},
                'Iowa':{'Illinois','Minnesota','Missouri','Nebraska','South_Dakota','Wisconsin'},
                'Kansas':{'Colorado','Missouri','Nebraska','Oklahoma'},
                'Kentucky':{'Illinois','Indiana','Missouri','Ohio','Tennessee','Virginia','West_Virginia'},
                'Louisiana':{'Arkansas','Mississippi','Texas'},
                'Maine':{'New_Hampshire'},
                'Maryland':{'Delaware','Pennsylvania','Virginia','West_Virginia'},
                'Massachusetts':{'Connecticut','New_Hampshire','New_York','Rhode_Island','Vermont'},
                'Michigan':{'Indiana','Ohio','Wisconsin'},
                'Minnesota':{'Iowa','North_Dakota','South_Dakota','Wisconsin'},
                'Mississippi':{'Alabama','Arkansas','Louisiana','Tennessee'},
                'Missouri':{'Illinois','Iowa','Arkansas','Kansas','Kentucky','Nebraska','Oklahoma','Tennessee'},
                'Montana':{'North_Dakota','South_Dakota','Idaho','Wyoming'},
                'Nebraska':{'Colorado','Iowa','Kansas','Missouri','South_Dakota','Wyoming'},
                'Nevada':{'Arizona','California','Idaho','Oregon','Utah'},
                'New_Hampshire':{'Massachusetts','Vermont','Maine'},
                'New_Jersey':{'Delaware','New_York','Pennsylvania'},
                'New_Mexico':{'Arizona','Colorado','Oklahoma','Texas','Utah'},
                'New_York':{'Connecticut','Massachusetts','New_Jersey','Pennsylvania','Vermont'},
                'North_Carolina':{'Georgia','South_Carolina','Tennessee','Virginia'},
                'North_Dakota':{'Minnesota','Montana','South_Dakota'},
                'Ohio':{'Indiana','Kentucky','Michigan','Pennsylvania','West_Virginia'},
                'Oklahoma':{'Arkansas','Colorado','Kansas','Missouri','New_Mexico','Texas'},
                'Oregon':{'California','Idaho','Nevada','Washington'},
                'Pennsylvania':{'Delaware','Maryland','New_Jersey','New_York','Ohio','West_Virginia'},
                'Rhode_Island':{'Connecticut','Massachusetts'},
                'South_Carolina':{'Georgia','North_Carolina'},
                'South_Dakota':{'Iowa','Minnesota','Montana','Nebraska','North_Dakota','Wyoming'},
                'Tennessee':{'Alabama','Arkansas','Georgia','Kentucky','Mississippi','Missouri','North_Carolina','Virginia'},
                'Texas':{'Arkansas','Louisiana','New_Mexico','Oklahoma'},
                'Utah':{'Arizona','Colorado','Idaho','Nevada','New_Mexico','Wyoming'},
                'Vermont':{'Massachusetts','New_Hampshire','New_York'},
                'Virginia':{'Kentucky','Maryland','North_Carolina','Tennessee','West_Virginia'},
                'Washington':{'Idaho','Oregon'},
                'West_Virginia':{'Kentucky','Maryland','Ohio','Pennsylvania','Virginia'},
                'Wisconsin':{'Illinois','Iowa','Michigan','Minnesota'},
                'Wyoming':{'Colorado','Idaho','Montana','Nebraska','South_Dakota','Utah'},  
    } 

    @staticmethod
    def complete(number_of_nodes):
        return {x : {y for y in range(number_of_nodes) if y != x} for x in range(number_of_nodes)}

    g1 = {  
            1: {2,3},
            2: {1,3},
            3: {1,2,4},
            4: {3,5,6},
            5: {4,6},
            6: {4,5,7},
            7: {6},
            }


if __name__ == '__main__':
    #print(america)

    graph = Graph(complete)
    print(graph.chromatic_number)
    print([graph.chromatic_polynomial(x) for x in range(10)])