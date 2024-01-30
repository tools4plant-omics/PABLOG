#Class to define Region object:
#in this application a Region is a segment of nucleotide bases, that starts at index s and ends at index e ([s,e])
#PABLOG tool identifies some Regions as potential primers sites (if found)  

class Region:
    
   def __init__(self, s, e, name="", cnt=1, cnf=0.0): 
      self.refname = name
      self.start = s
      self.end = e
      self.counter = cnt                           
      self.confidence = cnf                         

   def increase_counter(self):                     
      self.counter += 1

   def set_confidence(self, value):
      self.confidence = value
   
   def get_counter(self):  
      return str(self.counter)                     
   
   def __eq__(self, other):
      if (self.start == other.start) and (self.end == other.end):
         return True
      return False
   
   def __str__(self):
      string = '{}\t{}\t{}\t{}'.format(self.refname.ljust(10),str(self.start).ljust(10),str(self.end).ljust(10),str(round(self.confidence,2)).ljust(10))
      return string

def addUniqueRegion(unique_regions, region):

   if region not in unique_regions:
      unique_regions.append(region)
   else:
      index = unique_regions.index(region)
      unique_regions[index].increase_counter()
   

   
