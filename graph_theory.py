class Graph:

    chromatic_polynomial = None

    def __init__(self, adjacency, cache=None,level=1):
        self.adjacency = dict(adjacency)
        self.level = level
        self.cache = cache

    @property
    def vertices(self):
        return set(self.adjacency.keys())
    
    @property
    def edges(self):
        edges = 0
        for x in self.adjacency.values():
            edges += len(x)
        return edges
    
    @property
    def is_tree(self):
        if (self.edges/2) == len(self.vertices) - 1:
            return True
        else:
            return False
    @property
    def is_connected(self):
        vertex = tuple(self.vertices)[0]
        if self.traverse({vertex},vertex) == self.vertices:
            return True
        else:
            return False
    
    @staticmethod
    def subgraph_of_common_vertex(graph, vertex):
        v_subset = graph.adjacency[vertex].union({vertex})
        
        return Graph({v : {x for x in graph.adjacency[v] if x in v_subset} for v in v_subset})

    @staticmethod
    def subgraph(graph,vertices):
        return Graph({v : {x for x in graph.adjacency[v] if x in vertices} for v in vertices})

    def traverse(self, vertex_set, vertex):
        connected = self.adjacency[vertex]
        to_be_traversed = set([v for v in connected if v not in vertex_set])
        
        if len(to_be_traversed) == 0:
            return vertex_set
        else:
            for v in to_be_traversed:
                vertex_set = self.traverse(vertex_set.union(to_be_traversed),v)
            return vertex_set
        
    def del_contract_nodes(self):

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

                    
    def deletion_graph(self,vertex1,vertex2):
        
        adjacency = dict(self.adjacency)
        adjacency[vertex2] = set([v for v in adjacency[vertex2] if v != vertex1])
        adjacency[vertex1] = set([v for v in adjacency[vertex1] if v != vertex2])

        return Graph(adjacency=adjacency,level=self.level + 1)
        # return Graph(adjacency=adjacency,level=self.level + 1,cache=self.cache)
    
    def contraction_graph(self,vertex1,vertex2):

        adjacency = dict(self.adjacency)
        adjacency[vertex2] = set([v for v in adjacency[vertex2] if v != vertex1])
        adjacency[vertex1] = set([v for v in adjacency[vertex1] if v != vertex2]).union(adjacency[vertex2])

        for node in adjacency[vertex2]:
            adjacency[node] = set([v for v in adjacency[node] if v != vertex2]).union({vertex1})
            
        adjacency[vertex1].union(adjacency[vertex2])
        del adjacency[vertex2]

        return Graph(adjacency=adjacency,level=self.level + 1)
        # return Graph(adjacency=adjacency,level=self.level + 1,cache=self.cache)
    
    @property
    def deletion_contraction(self):
        
        if self.is_connected == False:
            print('does not exist: graph is not connected')
            return None
    
        if self.is_tree:
            
            return lambda n: n*(n-1)**(len(self.vertices) - 1)
        
        
        else:

            vertex1, vertex2 = self.del_contract_nodes()
            #deletion_poly    = self.cache.search_cache(self.deletion_graph(vertex1, vertex2))
            #contraction_poly = self.cache.search_cache(self.contraction_graph(vertex1, vertex2))
            deletion_poly    = self.deletion_graph(vertex1, vertex2).deletion_contraction
            contraction_poly = self.contraction_graph(vertex1, vertex2).deletion_contraction      
            return lambda G: deletion_poly(G) - contraction_poly(G)
            #return lambda G: self.cache.search_cache(self.deletion_graph(vertex1, vertex2))(G) - self.cache.search_cache(self.contraction_graph(vertex1, vertex2))(G)
    
    def glue1(self):

        vertices_remaining = [v for v in self.vertices]
        vertex_1 = vertices_remaining.pop()

        graph_main = self.subgraph_of_common_vertex(self,vertex_1)

        previous_main_chromatic_polynomial = graph_main.deletion_contraction



        while vertices_remaining:

            vertex_2 = [v for v in graph_main.vertices if v in vertices_remaining][0]
            vertices_remaining = [v for v in vertices_remaining if v != vertex_2]

            subgraph_vertex_2       = self.subgraph_of_common_vertex(self,vertex_2)
            intersection_subgraph   = self.subgraph(self,graph_main.vertices.intersection(subgraph_vertex_2.vertices))
            
            subgraph_vertex_2.chromatic_polynomial      = subgraph_vertex_2.deletion_contraction
            intersection_subgraph.chromatic_polynomial  = intersection_subgraph.deletion_contraction
            
            graph_main = self.subgraph(self,graph_main.vertices.union(subgraph_vertex_2.vertices))

            chromatic_polynomial = lambda G:  previous_main_chromatic_polynomial(G)*subgraph_vertex_2.chromatic_polynomial(G)/max(intersection_subgraph.chromatic_polynomial(G),1)
        


        self.chromatic_polynomial = chromatic_polynomial

           # previous_main_chromatic_polynomial = new_poly

        return [self.chromatic_polynomial(x) for x in range(10)]

    def glue(self):
        vertices_remaining = [v for v in self.vertices]
        vertex_1 = vertices_remaining.pop()

        initial_graph = self.subgraph_of_common_vertex(self,vertex_1)
        initial_graph.chromatic_polynomial = initial_graph.deletion_contraction

        self.chromatic_polynomial = lambda G: initial_graph.deletion_contraction(G)*self._glue(initial_graph, vertices_remaining)(G)



    def _glue(self, graph_main, vertices_remaining):

        print('here')

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

    @property
    def chromatic_number(self):
        if self.chromatic_polynomial == None:
            self.glue()

        chromatic_number = 1

        while self.chromatic_polynomial(chromatic_number) <= 0:
            chromatic_number += 1
        return chromatic_number

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

complete = {x : {y for y in range(8) if y != x} for x in range(8)}


if __name__ == '__main__':
    #print(america)
    graph = Graph(america)
    print(graph.chromatic_number)