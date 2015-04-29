# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import patches, cm, pyplot

class Model(object):
  """
  Problème de treillis
  """
  def __init__(self):
    self.nodes = []
    self.bars = []
    #self.blocked_dof = []
    #self.forces = []
  
  def add_node(self, node, *args, **kwargs):
    if isinstance(node, Node) == False:
      node = Node(node, *args, **kwargs)
    self.nodes.append(node)  
    return node
  
  def add_bar(self,*args , **kwargs):
    bar = Bar(*args , **kwargs)
    self.bars.append(bar)
    return bar
  
  def __repr__(self):
    return "<Model: {0} nodes, {1} bars>".format(len(self.nodes), len(self.bars))
    
  def get_stiffness_matrix(self):
    nodes = np.array(self.nodes)
    bars = self.bars
    nn = len(nodes)
    nb = len(bars)
    K = np.zeros([2 * nn, 2 * nn])
    for b in bars:
      conn = b.conn
      i0 = np.where(nodes == conn[0])[0][0]
      i1 = np.where(nodes == conn[1])[0][0]
      Kb = b.stiffness_matrix
      Kb_bar = Kb[0:2,0:2]
      K[2*i0:2*i0 +2,2*i0:2*i0 +2] += Kb_bar
      K[2*i1:2*i1 +2,2*i1:2*i1 +2] += Kb_bar
      K[2*i0:2*i0 +2,2*i1:2*i1 +2] -= Kb_bar
      K[2*i1:2*i1 +2,2*i0:2*i0 +2] -= Kb_bar
    return K
    
  stiffness_matrix = property(get_stiffness_matrix)   
  
  
  
  def add_force(self, node, magnitude):
    magnitude = np.array(magnitude).astype(np.float64)[0:2]
    self.forces.append( (node, magnitude) )
  
  def get_force_vector(self):
    nodes = self.nodes
    force_vector = np.array([n.force for n in nodes]).flatten()
    return force_vector
  force_vector = property(get_force_vector)
  
  def solve(self):
    adof = self.active_dof
    nodes= self.nodes
    nn = len(nodes)
    u = np.zeros(2 * nn)
    K = self.stiffness_matrix
    Adof1, Adof0 = np.meshgrid(adof, adof)
    Kr = K[(Adof0, Adof1)]
    f = self.force_vector
    fr = f[adof]
    u[adof] = np.linalg.solve(Kr,fr)
    nodes = np.array(self.nodes)
    f = np.dot(K, u)
    for i in xrange(len(nodes)):
      node = nodes[i]
      for j in xrange(2):
        node.displacement[j] = u[2*i+j]
        if node.block[j]: node.force[j] = f[2*i+j]
          
      
    bars = self.bars
    for bar in bars:
      n0 = bar.conn[0]
      n1 = bar.conn[1]
      bar.elongation = (n1.displacement - n0.displacement).dot(bar.direction())  
      bar.tension = bar.elongation * bar.stiffness
      bar.strain = bar.elongation / bar.length()
      bar.stress = bar.tension / bar.section 
   
  def get_blocked_dof(self):
    bdof = []
    nodes = self.nodes
    for node in nodes:
      if node.block[0]: bdof.append((node, 0))
      if node.block[1]: bdof.append((node, 1))   
    return bdof
  blocked_dof = property(get_blocked_dof)
  
  def get_active_dof(self):
    nodes = self.nodes
    nn = len(nodes)
    a = np.ones(2 * nn)
    for i in xrange(nn):
      if nodes[i].block[0]: a[2 * i    ] = 0 
      if nodes[i].block[1]: a[2 * i + 1] = 0  
    a = np.where(a == 1)[0]
    return a
  active_dof = property(get_active_dof)    
  
  def bbox(self, deformed = True, factor = .2):
    xlim = np.zeros(2)
    ylim = np.zeros(2)
    for n in self.nodes:
      pos = n.coords.copy()
      if deformed: pos += n.displacement
      xlim[0] = min(xlim[0], pos[0])
      xlim[1] = max(xlim[1], pos[0])
      ylim[0] = min(ylim[0], pos[1])
      ylim[1] = max(ylim[1], pos[1])
    d = max(xlim[1]-xlim[0], ylim[1]-ylim[0])   
    xlim[0] -= d*factor
    xlim[1] += d*factor
    ylim[0] -= d*factor
    ylim[1] += d*factor
    return xlim, ylim
  
  def draw(self, ax, deformed = True, field = "stress", label = True, forces = True, displacements = False, force_scale = 1., displacement_scale = 1.):
    for node in self.nodes: node.draw(ax, deformed = deformed, label = label)
    bars = self.bars
    length = np.array([b.length for b in bars])
    if field != None:
      if field == "stress": 
        values = np.array([b.stress for b in bars])
        cbar_label = "Normal stress, $\sigma$"
      if field == "tension": 
        values = np.array([b.tension for b in bars])
        cbar_label = "Tension, $N$"
      colors = cm.jet(values)
      vmin = min(0., values.min())
      vmax = max(0., values.max())
      colormap = pyplot.cm.ScalarMappable(cmap=cm.jet, 
        norm=pyplot.Normalize(vmin = vmin, vmax = vmax))
      colormap._A = []
    for i in xrange(len(bars)):
      if field == None:
        color = None
      else:  
        color = colormap.to_rgba(getattr(bars[i], field))
      bars[i].draw(ax = ax, deformed = deformed, color = color)
    if field != None:
      cbar = pyplot.colorbar(colormap)
      cbar.set_label(cbar_label)    
    F = np.array([node.force for node in self.nodes]).transpose()
    P = np.array([node.coords for node in self.nodes]).transpose()
    U = np.array([node.displacement for node in self.nodes]).transpose()
    if deformed : P += U
    if forces:
      qf = ax.quiver(P[0], P[1], F[0], F[1], scale_units='xy', angles = "xy", pivot="tail", scale=force_scale**-1, color = "red")
      #qk = ax.quiverkey(qf, 0.1, 1.1, 1, r'1 N', labelpos='E')
    if displacements:
      if deformed: 
        upos = "tip"
      else:
        upos = "tail"  
      qu = ax.quiver(P[0], P[1], U[0], U[1], scale_units='xy', angles = "xy", pivot=upos, scale=1., color = "green")
  
  def get_mass(self):
    return np.array([b.mass for b in self.bars]).sum()
  mass = property(get_mass)  
      
class Node(object):
  def __init__(self, coords = np.array([0., 0.]), label = None, bock_side = 1):
    coords = np.array(coords).astype(np.float64)[0:2]
    self.coords = np.array(coords)
    self.displacement = np.zeros(2)
    self.force = np.zeros(2)
    self.block = np.array([False, False])
    self.label = label
    self.block_side = block_side
  
  def __repr__(self):
    return "<Node {0}: x = {1}, y = {2}>".format(self.label, self.coords[0], self.coords[1])
  
  def draw(self, ax, deformed = True, radius = 0.03, label = True, force_factor = 5.):
    pos = self.coords.copy()
    if deformed: pos += self.displacement
    patch = patches.Circle(pos, radius, color='k',clip_on=False)
    ax.add_artist(patch)
    if label == True:
      if self.label != None:
        an = ax.annotate(self.label, xy=pos,  xycoords='data',
                    xytext=(-50, 30), textcoords='offset points',
                    bbox=dict(boxstyle="round", fc="1.",clip_on=False),
                    arrowprops=dict(arrowstyle="<-",
                                    connectionstyle="angle,angleA=0,angleB=90,rad=10")
                    )
        an.draggable()
    if self.block[0]:
      d = radius * 2.
      bs = self.block_side
      verts = np.array([[-.1,0.], [-1.,.9], [-1.,-.9], [-.1, 0.]]) 
      verts *= d
      verts += pos
      p = patches.Polygon(verts, facecolor = "none",clip_on=False, linewidth = 1.5)
      ax.add_artist(p)
      p = patches.Circle(pos + np.array([-1.5, .5])*d, d*.5, facecolor='none',clip_on=False, linewidth = 1.5)  
      ax.add_artist(p) 
      p = patches.Circle(pos + np.array([-1.5, -.5])*d, d*.5, facecolor='none',clip_on=False, linewidth = 1.5)  
      ax.add_artist(p) 
    if self.block[1]:
      d = radius * 2.
      
      verts = np.array([[0.,-.1], [-.9, -1.], [.9,-1.], [0., -.1]]) 
      verts *= d
      verts += pos
      p = patches.Polygon(verts, facecolor = "none",clip_on=False, linewidth = 1.5)
      ax.add_artist(p)
      p = patches.Circle(pos + np.array([-.5, -1.5])*d, d*.5, facecolor='none',clip_on=False, linewidth = 1.5)  
      ax.add_artist(p) 
      p = patches.Circle(pos + np.array([.5, -1.5])*d, d*.5, facecolor='none',clip_on=False, linewidth = 1.5)  
      ax.add_artist(p)   
    
      
      
    
class Bar(object):
  def __init__(self, n1, n2, section = 1., modulus = 1., density = 1.):
    self.conn = [n1, n2]
    self.section = float(section)
    self.modulus = float(modulus)
    self.density = float(density)
    self.tension = 0.
    self.elongation = 0.
    self.strain = 0.
    self.stress = 0.
    
  def __repr__(self):
    return "<Bar: ({0}, {1}) -> ({2}, {3}), S = {4}, E = {5}, rho = {6}>".format(
      self.conn[0].coords[0], 
      self.conn[0].coords[1],
      self.conn[1].coords[0],
      self.conn[1].coords[1], 
      self.section, self.modulus, self.density)  
  
  def draw(self, ax, deformed = True, offset = .1, width = .05, color = None):
    b = self
    o = offset
    w = width
    l = b.length(deformed)
    n0,n1 = b.conn[0].coords.copy(), b.conn[1].coords.copy()
    if deformed:
      n0 += b.conn[0].displacement
      n1 += b.conn[1].displacement
    u = b.direction(deformed)
    n = b.normal(deformed)
    verts = np.array([ n0, n0 + o*u, n0 + (o+w)*u + w*n, 
      n0 + (l-o-w)*u + w*n, n0 + (l-o)*u, n0 + l*u, n0 + (l-o)*u,
      n0 + (l-o-w)*u - w*n, n0 + (o+w)*u - w*n, n0 + o*u, n0])
    tension = b.tension
    if color == None: color = "none"
    p = patches.Polygon(verts, facecolor = color, linewidth = 1.5,clip_on=False)
    ax.add_artist(p)  
   
  def length(self, deformed = False):
    conn = self.conn
    if deformed:
      return ((conn[0].coords + conn[0].displacement 
        - conn[1].coords - conn[1].displacement)**2).sum()**.5 
    else:
      return ((conn[0].coords - conn[1].coords)**2).sum()**.5 
   
     
  def get_volume(self):
    S = self.section
    L = self.length()
    return S * L
  volume = property(get_volume)   
  
  def get_mass(self):
    V = self.volume
    rho = self.density
    return V * rho
  mass = property(get_mass)   
  
  def get_stiffness(self):
    S = self.section
    L = self.length()
    E = self.modulus
    return E * S / L    
  stiffness = property(get_stiffness)   
  
  def direction(self, deformed = False):
    conn = self.conn
    if deformed:
      return (conn[1].coords + conn[1].displacement 
        - conn[0].coords - conn[0].displacement) / self.length(deformed)
    else:
      return (conn[1].coords - conn[0].coords) / self.length(deformed)
  
  
  def normal(self, deformed = False):
    t = self.direction(deformed)
    n = np.zeros(2)
    n[0] = -t[1]
    n[1] = t[0]
    return n
  
  
  
  def get_stiffness_matrix(self):
    u = self.direction()
    k = self.stiffness
    ux, uy = u[0], u[1]
    a = np.array([[-ux, -uy, ux, uy]])
    K = k *  a.transpose().dot(a)
    return K
  stiffness_matrix = property(get_stiffness_matrix)   
          
