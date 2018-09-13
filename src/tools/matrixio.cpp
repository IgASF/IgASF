/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <fcntl.h>

#include <tools/matrixio.h>

outData::outData(const RowM &mat )
: rowMajor(true), rows(mat.rows()), cols(mat.cols()), nnzs(mat.nonZeros()),
innerBeg(mat.outerIndexPtr()),innerPos(mat.innerIndexPtr()),values(mat.valuePtr())
{}

outData::outData(RMap mat )
: rowMajor(true), rows(mat.rows()), cols(mat.cols()), nnzs(mat.nonZeros()),
innerBeg(mat.outerIndexPtr()),innerPos(mat.innerIndexPtr()),values(mat.valuePtr())
{}

outData::outData(const ColM &mat )
: rowMajor(false), rows(mat.rows()), cols(mat.cols()), nnzs(mat.nonZeros()),
innerBeg(mat.outerIndexPtr()),innerPos(mat.innerIndexPtr()),values(mat.valuePtr())
{}

outData::outData(CMap mat )
: rowMajor(false), rows(mat.rows()), cols(mat.cols()), nnzs(mat.nonZeros()),
innerBeg(mat.outerIndexPtr()),innerPos(mat.innerIndexPtr()),values(mat.valuePtr())
{}

RMap inData::asRowMatrix() const
{
    if (rowMajor) return RMap(rows,cols,values.size(), innerBeg.data(), innerPos.data(), values.data());
    else return RMap(0,0,0, nullptr,nullptr,nullptr);
}
CMap inData::asColMatrix() const
{
    if (rowMajor) return CMap(0,0,0, nullptr,nullptr,nullptr);
    else return CMap(rows,cols,values.size(), innerBeg.data(), innerPos.data(), values.data());
}
bool  inData::valid() const
{
    return rows>0 && cols>0;
}
double  inData::norm() const
{
    return rowMajor ?  asRowMatrix().norm() : asColMatrix().norm();
}

namespace {

struct resultHeader
{
    int    rowMajor;
    int    rows;
    int    cols;
    int    nnzs;
};

int fd_set_blocking(int fd, int blocking) {
    /* Save the current flags */
    int flags = fcntl(fd, F_GETFL, 0);
    if (flags == -1)
        return 0;
    
    if (blocking)
        flags &= ~O_NONBLOCK;
    else
        flags |= O_NONBLOCK;
    return fcntl(fd, F_SETFL, flags) != -1;
}

void writeN (const char *start,unsigned long limit,int fd)
{
    long cur;
    unsigned long written=0;
    while ( written < limit )
        {
        unsigned long clipped=std::min<unsigned long>(limit-written,  INT_MAX );
        cur=write(fd,start+written, clipped); // the limit is because of a MACOS bug in POSIX
        if (cur<0)
            {
            int err=errno;
            switch (err) {
                case ERANGE:      throw std::runtime_error("Error writing: chunk is too big (ERANGE)");
                case EPIPE:       throw std::runtime_error("Error writing: pipe/socket not connected for writing (EPIPE)");
                case ENOSPC:      throw std::runtime_error("Error writing: not enough space (ENOSPC)");
                case EFBIG:       throw std::runtime_error("Error writing: file is too big (EFBIG)");
                case EBADF:       throw std::runtime_error("Error writing: invalid fs (EBADF)");
                case EIO:         throw std::runtime_error("Error writing: physical-fail or write to controlling terminal (EIO)");
                case ECONNRESET:  throw std::runtime_error("Error writing: socket not connected (ECONNRESET)");
                case EINVAL:      throw std::runtime_error("Error writing: negative offset or fd linked to multiplexer (EINVAL)");
                case EFAULT:      throw std::runtime_error("Error writing: illegal address (EFAULT)");
                case ENOBUFS:     throw std::runtime_error("Error writing: not enough system resources (ENOBUFS)");
                case ENXIO:       throw std::runtime_error("Error writing: nonexistent device (ENXIO)");
                case EACCES:      throw std::runtime_error("Error writing: no write permission (EACCES)");
                case ENETDOWN:    throw std::runtime_error("Error writing: network is down (ENETDOWN)");
                case ENETUNREACH: throw std::runtime_error("Error writing: network destination unreachable (ENETUNREACH)");
                case EINTR:       continue;
                case EAGAIN: fd_set_blocking(fd, true); break; // force blocking IO
                default: throw std::runtime_error("Error writing: other error");
            }
            }
        written +=cur;
        }
}

void readN( char *start, unsigned long limit, int fd )
{
    long cur;
    unsigned long done = 0;
    while ( done < limit )
        {
        unsigned long clipped=std::min<unsigned long>(limit-done,  INT_MAX ); // the limit is because of a MACOS bug in POSIX
        cur=read(fd,start+done, clipped);
        if      (cur<0 && errno==EINTR) continue;
        else if (cur<0 && errno!=EINTR) throw std::runtime_error("Error while reading: errno = "+std::to_string(errno));
        if (cur==0)                     throw std::runtime_error("Error while reading: unexpected end of file.");
        done += cur;
        }
}

} // anonymous namespace

void writeMatrix(const outData &data, const std::string& fn)
{
    int fd;
    if (fn == "stdout")
        fd = STDOUT_FILENO;
    else
        fd = open(fn.c_str(),O_CREAT|O_TRUNC|O_WRONLY,0644);
    
    size_t limit = 0;
    long matrixSize = data.rowMajor ? data.rows+1 : data.cols+1;
    long nnzs = data.nnzs;
    
    // write lengths in a header
    resultHeader head={
        data.rowMajor,
        data.rows,
        data.cols,
        data.nnzs
    };
    const char *start=reinterpret_cast<char*>(&head);
    writeN(start,sizeof(head),fd);
    
    // write null byte
    char nullByte = '\0';
    writeN(&nullByte, 1, fd);

    // write begin of inner vectors
    limit = sizeof(int)*matrixSize;
    start = reinterpret_cast<const char*>(data.innerBeg);
    writeN(start,limit,fd);
    
    // write position of coefficients
    limit = sizeof(int)*nnzs;
    start = reinterpret_cast<const char*>(data.innerPos);
    writeN(start,limit,fd);
    
    // write coefficients
    limit = sizeof(double)*nnzs;
    start = reinterpret_cast<const char*>(data.values);
    writeN(start,limit,fd);
    
    if (fn != "stdout")
        close(fd);
}

inData readMatrix(const std::string& fn)
{
    int fd;
    if (fn == "stdin")
        fd = STDIN_FILENO;
    else
        fd = open(fn.c_str(),O_RDONLY);
    
    size_t limit;
    resultHeader head;
    inData      result;
    
    // read header with lengths
    char *start = reinterpret_cast<char*>(&head);
    readN(start, sizeof(head), fd);
    
    result.rowMajor=head.rowMajor;
    result.cols=head.cols;
    result.rows=head.rows;
    
    // read null byte
    char nullByte = 'a';
    readN(&nullByte, 1, fd);
    if (nullByte != '\0') throw std::runtime_error("No nullbyte found.");

    result.innerBeg.resize( head.rowMajor ? head.rows+1 : head.cols +1 );
    result.innerPos.resize( head.nnzs );
    result.values.resize( head.nnzs );

    // read begin of inner vectors
    limit = sizeof(int)*result.innerBeg.size();
    start = reinterpret_cast<char*>(result.innerBeg.data());
    readN(start, limit, fd);
    
    // read position of coefficients
    limit = sizeof(int)*result.innerPos.size();
    start = reinterpret_cast<char*>(result.innerPos.data());
    readN(start, limit, fd);
    
    // read coefficients
    limit = sizeof(double)*result.values.size();
    start = reinterpret_cast<char*>(result.values.data());
    readN(start, limit, fd);
    
    if (fn != "stdin")
        close(fd);
    
    return result;
}
