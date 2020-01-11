##TGC - **T**housands **G**enome **C**ompressor

###Getting Started 

```
# Clone repository
git clone https://github.com/refresh-bio/tgc.git
# Build 
cd tgc
cd src
make
cd ..
# Repeat experiments described in the publication
./run -abcdefgijklcd tgc
```


### TGC—What is it?

Thousands Genome Compressor is a tool to estimate the boundaries of compression ratio for human genome compression. It can be also used as a very effective tool for compression Variant Call Format (VCF) files.

#### The architecture of TGC

TGC is composed of several programs that were used in our experiments on the genomes from the 1000 Genomes project. A description of the tool can be found in TGC description.

#### How good is TGC?

We were able to compress the genomes from the 1000GP about 15,500 times. More details can be found in our paper pointed below.

### Terms of use of TGC

TGC is in general a free compression program available in source code release. 

### Publications

 + Deorowicz, S., Danek, A., Grabowski, Sz., Genome compression: a novel approach for large collections, [Bioinformatics](http://bioinformatics.oxfordjournals.org/content/29/20/2572), 2013; 29(20):2572–2578.


### Developers

The TGC algorithm was invented by Sebastian Deorowicz, Szymon Grabowski, and Agnieszka Danek.
The implementation is by [Agnieszka Danek](https://github.com/agnieszkadanek) and Sebastian Deorowicz.
