FROM kbase/kbase:sdkbase.latest
MAINTAINER Dylan Chivian
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.



# Update Transform (should go away eventually)
RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf transform && \ 
  git clone https://github.com/kbase/transform -b develop

# setup the transform, but ignore errors because sample data cannot be found!
RUN \
  . /kb/dev_container/user-env.sh; \
  cd /kb/dev_container/modules/transform/t/demo; \
  python setup.py; \
  exit 0;



# RUN apt-get update
# -----------------------------------------

# Install SDK Module
#
RUN mkdir -p /kb/module
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
WORKDIR /kb/module
RUN make


# Install VSEARCH
#
RUN mkdir -p /kb/module/vsearch
WORKDIR /kb/module
RUN git clone https://github.com/torognes/vsearch
WORKDIR vsearch
RUN ./configure 
RUN make
RUN make install


# add test data
#
RUN mkdir -p /kb/module/test
WORKDIR /kb/module/test
RUN git clone https://github.com/dcchivian/kb_vsearch_test_data
WORKDIR /kb/module


ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]


