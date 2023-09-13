FROM ubuntu:20.04
RUN apt-get update -y
RUN apt-get install python3-dev -y
RUN apt-get install python3-pip -y
RUN apt-get install python3-venv -y
RUN apt-get install git -y
RUN apt-get install libhdf5-dev -y
RUN apt-get install curl -y
RUN apt-get install libbz2-dev -y
RUN apt-get install liblzma-dev -y
RUN apt-get autoremove -y

# --- prep for hdf5 for HelixerPost --- #
WORKDIR /tmp/
RUN curl -L https://github.com/h5py/h5py/releases/download/3.2.1/h5py-3.2.1.tar.gz --output h5py-3.2.1.tar.gz
RUN tar -xzvf h5py-3.2.1.tar.gz
#RUN gcc -O2 -fPIC -shared -Ilzf -I/usr/include/hdf5/serial/ lzf/*.c lzf_filter.c -lhdf5 -L/lib/x86_64-linux-gnu/hdf5/serial -o liblzf_filter.so
RUN cd h5py-3.2.1/lzf/ && gcc -O2 -fPIC -shared -Ilzf -I/usr/include/hdf5/serial/ lzf/*.c lzf_filter.c  -lhdf5_cpp -L/lib/x86_64-linux-gnu/hdf5/serial -o liblzf_filter.so
RUN mkdir /usr/lib/x86_64-linux-gnu/hdf5/plugins && mv h5py-3.2.1/lzf/liblzf_filter.so /usr/lib/x86_64-linux-gnu/hdf5/plugins
RUN rm -r /tmp/h5py*

# --- rust install for HelixerPost --- # 
RUN curl https://sh.rustup.rs -sSf > rustup.sh
RUN chmod 755 rustup.sh
RUN ./rustup.sh -y
RUN rm rustup.sh
ENV PATH="/root/.cargo/bin:${PATH}"

# --- Helixer and HelixerPost --- #

WORKDIR /opt/helixer_user/
RUN git clone -b v0.3.0 https://github.com/weberlab-hhu/Helixer.git Helixer
RUN pip install --no-cache-dir -r /opt/helixer_user/Helixer/requirements.txt
RUN cd Helixer && pip install --no-cache-dir .

WORKDIR /opt/helixer_user/
RUN git clone https://github.com/TonyBolger/HelixerPost.git
RUN cd HelixerPost && git checkout f7cb814649150f01a17be2ee99fd57b027afd82d
RUN mkdir bin
ENV PATH="/opt/helixer_user/bin:${PATH}"
RUN cd HelixerPost/helixer_post_bin && cargo build --release
RUN mv /opt/helixer_user/HelixerPost/target/release/helixer_post_bin /opt/helixer_user/bin/
RUN rm -r /opt/helixer_user/HelixerPost/target/release/

RUN Helixer/scripts/fetch_helixer_models.py
#RUN sed -i "s/MODEL_PATH,/\"/root/.local/share/Helixer/models/\",/g" /usr/local/bin/Helixer.py  
# Install Scripts
ENV repodir /seqprg
RUN mkdir -p /seqprg;
ADD scripts /seqprg/scripts


ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
