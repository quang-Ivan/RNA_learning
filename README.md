# RNA_learning

Learning RNA and AI

# 构建docker 镜像的方法

```
docker build -t rnabio_enhanced .
```

# 运行的方法

```
docker run --user ubuntu:ubuntu -it rnabio_enhanced
```

```
docker run -it -v /E/RNA-data:/home/ubuntu/workspace/DATA aierlma/rnabio_enhanced:latest /bin/bash
```

底下是挂载路径
