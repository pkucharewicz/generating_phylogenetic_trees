import datetime
from typing import List, Tuple
from Bio import SeqIO
import datetime
from collections import defaultdict
from random import randint
import math

class Sample:
    sample_id: str
    country: str
    date: datetime.date
    seq: str    

    def __init__(self,sample_id: str,country: str,date: datetime.date,seq: str):
        self.sample_id = sample_id
        self.country = country
        self.date = date
        self.seq=seq

class Node:
    value: Sample
    children: List

    def __init__(self, value,children=None,parent=None):
        self.value = value
        if children==None:
            children=[]
        self.children = children
        self.parent=parent

    def add_child(self,node):
        self.children.append(node)

    def find_parent(self,node):
        min=edit_distance(self.value.seq,node.value.seq)
        parent=self
        if self.children!=None:
            for child in self.children:
                res=child.find_parent(node)
                if res[0]<min:
                    min=res[0]
                    parent=res[1]
        return(min,parent)


    def edges_node(self,list):
        if self.children!=None:
            for child in self.children:
                list.append((self.value.sample_id,child.value.sample_id))
                child.edges_node(list)

    def cost_node(self,n):
        if self.children!=None:
            for child in self.children:
                n+=edit_distance(child.value.seq,self.value.seq)
                n=child.cost_node(n)
        return n

    def filter_node(self,country,parent=None):
        res=[]
        if self.value.country==country:
            node=Node(self.value)
            if parent:
                parent.add_child(node)
                node.parent=parent
            else:
                res.append(Tree(node))
            parent=node

        if self.children!=None:
            for child in self.children:
                res.extend(child.filter_node(country,parent))

        return(res)

class Tree:

    def __init__(self, node=None):
        self.root = node

    def add_node(self,node):
        if self.root!=None:
            parent=self.root.find_parent(node)[1]
            parent.add_child(node)
            node.parent=parent
        else:
            self.root=node

    def edges(self) -> List[Tuple[str, str]]:
        '''returns a list of all edges in the tree in O(n) time'''
        pairs=[]
        if self.root!=None:
            self.root.edges_node(pairs)
        return(pairs)

    def filter(self, country: str) -> List['Tree']:
        '''returns a list of new trees that consist of samples from a given country in a O(n) time'''
        res=[]
        if self.root!=None:
            res=self.root.filter_node(country)
        return(res)

    def cost(self):
        res=0
        if self.root!=None:
            res=self.root.cost_node(res)

        return(res)

def edit_distance(a,b):
    dp = defaultdict(int)

    for i in range(len(a)):
        dp[i,0]=i

    for i in range(len(b)):
        dp[0,i]=i

    for j in range(1,len(b)+1):
        for i in range(1,len(a)+1):
            if a[i-1]==b[j-1]:
                substitution_cost=0
            else:
                substitution_cost=1

            dp[i, j] = min(dp[i-1, j] + 1,
                            dp[i, j-1] + 1,
                            dp[i-1, j-1] + substitution_cost)

    return dp[len(a), len(b)]

def get_distance(node1,node2,distances):
    try:
        return(distances[(node1.value.sample_id,node2.value.sample_id)])
    except:
        distance=edit_distance(node1.value.seq,node2.value.seq)
        distances[(node1.value.sample_id,node2.value.sample_id)]=distances[(node2.value.sample_id,node1.value.sample_id)]=distance
        return(distance)


def read_data(filename: str) -> List[Sample]:
    '''
    reads info about samples from a given fasta file, creates Sample objects and sorts them with increasing date of sample's collection.
    '''
    records = SeqIO.parse(filename, "fasta")
    samples=[]
    for r in records:
        data=r.description.split("|")
        id=data[0]
        country=data[1]
        date=datetime.datetime.strptime(data[2], '%Y-%m-%d').date()
        seq=str(r.seq)
        samples.append(Sample(id,country,date,seq))
        samples.sort(key=lambda sample: sample.date)
    return samples

def construct_optimal_tree(samples: List[Sample]) -> Tree:
    '''
    constructs optimal phylogenetic tree for a given list of samples in O(n^2 * m^2), where n is a no. of samples and m is a maximal sequence length
    '''
    tree=Tree()
    for s in samples:
        child=Node(s)
        tree.add_node(child)

    return tree

def approximate_parent(node,random_node,distances):
    min_distance=get_distance(random_node,node,distances)
    parent=random_node
    if random_node.children!=None:
        for child in random_node.children:
            distance=get_distance(child,node,distances)
            if distance<min_distance:
                min_distance=distance
                parent=child
    if random_node.parent!=None:
        distance=get_distance(random_node.parent,node,distances)
        if distance<min_distance:
            min_distance=distance
            parent=random_node.parent

    if parent==random_node:
        return(parent,min_distance)
    else:
        return(approximate_parent(node,parent,distances))

def construct_approximate_tree(samples: List[Sample]) -> Tree:
    '''
    constructs approximate phylogenetic tree for given samples
    '''
    distances={}
    root=Node(samples[0])
    tree=Tree(root)
    nodes=[root]
    for i in range(1,len(samples)):
        node=Node(samples[i])
        parent=None
        min_distance=math.inf
        for j in range(3):
            random_node=nodes[randint(0,i-1)]
            res_node,distance=approximate_parent(node,random_node,distances)
            if distance<min_distance:
                min_distance=distance
                parent=res_node
        node.parent=parent
        parent.add_child(node)
        nodes.append(node)
    return(tree)
